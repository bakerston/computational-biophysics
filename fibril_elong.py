import Bio
import pandas as pd
from pandas import DataFrame, Series
import numpy as np
import math
from biopandas.pdb import PandasPdb
import csv



'''
Read base file
Find center axis
Calculate rotate angel
Extend fibril
output to new pdb
'''

def get_cir(cor_list):
    x1,y1,x2,y2,x3,y3=cor_list[0],cor_list[1],cor_list[2],cor_list[3],cor_list[4],cor_list[5]
    A=x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2
    B=(x1*x1+y1*y1)*(y3-y2)+(x2*x2+y2*y2)*(y1-y3)+(x3*x3+y3*y3)*(y2-y1)
    C=(x1*x1+y1*y1)*(x2-x3)+(x2*x2+y2*y2)*(x3-x1)+(x3*x3+y3*y3)*(x1-x2)
    D=(x1*x1+y1*y1)*(x3*y2-x2*y3)+(x2*x2+y2*y2)*(x1*y3-x3*y1)+(x3*x3+y3*y3)*(x2*y1-x1*y2)
    x_cir=-B/(2*A)
    y_cir=-C/(2*A)

    r_cir=math.sqrt((B**2+C**2-4*A*D)/(4*A**2))
    return [x_cir,y_cir,r_cir]

def degree(a,b,c):
    ba=a-b
    bc=c-b
    cosine_angle=np.dot(ba,bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
    angle=np.arccos(cosine_angle)
    return np.degrees(angle)

def rotate_point(P0,P,degree):
    #use -1 when chains rotate clockwise, use + when it rotates anti-closewise.
    rad = -math.radians(degree)
    x0, y0 = P0[0], P0[1]
    p_in = np.matrix([[P[0]], [P[1]], [1]])
    M = np.matrix([[math.cos(rad), -math.sin(rad), (1 - math.cos(rad)) * x0 + math.sin(rad) * y0],
                   [math.sin(rad), math.cos(rad), (1 - math.cos(rad)) * y0 - math.sin(rad) * x0], [0, 0, 1]])
    t = np.matmul(M, p_in)
    #    print(M)
    #    print(math.cos(3.14))
    return np.asarray(t).reshape(-1)[:2]

def unit_vector(p):
    return p/np.linalg.norm(p)



def rotate_3d_point(P,P1,P2,degree):
    rad = -math.radians(degree)
    #base1 upper [1290:], base2 lower [:1290]
    P=np.array(P)
    P1=np.array(P1)
    P2=np.array(P2)
    x,y,z = P[0], P[1], P[2]
    x1,y1,z1 = P1[0],P1[1],P1[2]
    x2,y2,z2 = P2[0],P2[1],P2[2]

    U=unit_vector(np.array(P1-P2))
    a, b, c = U[0], U[1], U[2]
    d = math.sqrt(b**2+c**2)

    T = np.matrix([[1,0,0,-x1],[0,1,0,-y1],[0,0,1,-z1],[0,0,0,1]])
    Rx = np.matrix([[1,0,0,0],[0,c/d,-b/d,0],[0,b/d,c/d,0],[0,0,0,1]])
    Ry = np.matrix([[d,0,-a,0],[0,1,0,0],[a,0,d,0],[0,0,0,1]])
    Rz = np.matrix([[math.cos(rad),-math.sin(rad),0,0],[math.sin(rad),math.cos(rad),0,0],[a,0,1,0],[0,0,0,1]])
    

    res = np.linalg.multi_dot([np.linalg.inv(T),np.linalg.inv(Rx),np.linalg.inv(Ry),Rz,Ry,Rx,T])
    ans= np.matmul(res,np.matrix([[x],[y],[z],[1]]))
    
    return np.asarray(ans).reshape(-1)[:3]



def write_to_pdb(input_frame,file_name):
    f_out = open(file_name, 'w')
    for index, row in input_frame.iterrows():
        f_out.write((
                                '%-6s%5d' + '  ' + '%-4s'  + '%3s' + ' ' + '%1s' + '%4d' + '    ' + 3 * '%-8.3f' + 2 * '%-6.2f' + '      ' + '%4s' + '%2s' + '  ' + '\n') % (
                    row['record_name'], row['atom_number'], row['atom_name'], row['residue_name'], row['chain_id'],
                    row['residue_number'], row['x_coord'], row['y_coord'], row['z_coord'], row['occupancy'],
                    row['b_factor'], row['segment_id'], row['element_symbol']))
    f_out.write('TER')
    return

def modi_layer(new_base_1,ang,pos_mean,z_increase,index):
    for i, trial in new_base_1.iterrows():
        tmp = rotate_point(pos_mean, [trial['x_coord'], trial['y_coord']], -index*ang)
        rounded = list((map(lambda x: round(x, 3), tmp)))
        new_base_1.loc[i, 'x_coord'], new_base_1.loc[i, 'y_coord'] = rounded[0], rounded[1]

        cur_char = new_base_1.loc[i, 'chain_id']
        new_base_1.loc[i, 'chain_id'] = chr(ord(cur_char) + 2*index)
    new_base_1['z_coord'] = new_base_1['z_coord'] + index*z_increase
    return

def translate_layer(new_base,trans_vec,index):
    new_base['x_coord'] = new_base['x_coord'] + trans_vec[0]*index
    new_base['y_coord'] = new_base['y_coord'] + trans_vec[1]*index
    new_base['z_coord'] = new_base['z_coord'] + trans_vec[2]*index
    return

def rotate_layer(new_base,ang,P1,P2,index):
    for i, trial in new_base.iterrows():
        tmp=rotate_3d_point([trial['x_coord'], trial['y_coord'],trial['z_coord']],P1,P2,ang*index)
        rounded = list((map(lambda x: round(x,3),tmp)))
        new_base.loc[i, 'x_coord'], new_base.loc[i, 'y_coord'] = rounded[0], rounded[1]
        cur_char = new_base.loc[i, 'chain_id']
        new_base.loc[i, 'chain_id'] = chr(ord(cur_char) + 2*index)
    return

'''
read pdb to dataframe
'''
base=PandasPdb().read_pdb('base.pdb')
starter=pd.DataFrame(base.df['ATOM'])
x=starter[1290:]


#base1 the upper layer
base1=starter[1290:]
chain3=starter[1290:1935]
chain4=starter[1935:]
#base2 the lower layer
base2=starter[:1290]
chain1=starter[:645]
chain2=starter[645:1290]

"""
chain 1~4

fibril1       fibril2
 -----3      ------4        :layer1  upper 
-----1        -----2        :layer2  lower
"""

"""
1.Find axis equation for fibril 1 and fibril 2.
2.compute translation of each fibril by one layer, delta_x,y,z
3.compute rotate angle of each fibril by one layer, theta.

"""

coord=DataFrame(starter,columns=['atom_name','x_coord','y_coord','z_coord'])
coord1=DataFrame(base1,columns=['atom_name','x_coord','y_coord','z_coord'])
coord2=DataFrame(base2,columns=['atom_name','x_coord','y_coord','z_coord'])

coord_c1=DataFrame(chain1,columns=['atom_name','x_coord','y_coord','z_coord'])
coord_c2=DataFrame(chain2,columns=['atom_name','x_coord','y_coord','z_coord'])
coord_c3=DataFrame(chain3,columns=['atom_name','x_coord','y_coord','z_coord'])
coord_c4=DataFrame(chain4,columns=['atom_name','x_coord','y_coord','z_coord'])



xy_mean=np.array(coord.mean(axis=0))

base1_mean=np.array(coord1.mean(axis=0))
base2_mean=np.array(coord2.mean(axis=0))

c1_mean=np.array(coord_c1.mean(axis=0))
c2_mean=np.array(coord_c2.mean(axis=0))
c3_mean=np.array(coord_c3.mean(axis=0))
c4_mean=np.array(coord_c4.mean(axis=0))

print(c3_mean)
print(c3_mean[:2])

"""
print(xy_mean)
print("base1 mean is: ", base1_mean)
print("base2 mean is:", base2_mean)
print("c1 mean is:", c1_mean)
print("c2 mean is:", c2_mean)
print("c3 mean is:", c3_mean)
print("c4 mean is:", c4_mean)
"""
"""
translation vector
"""
translation_1=c3_mean-c1_mean
translation_2=c4_mean-c2_mean
trans_vec_1=list(map(lambda x: round(x,3),translation_1))
trans_vec_2=list(map(lambda x: round(x,3),translation_2))

print(trans_vec_1)
print(trans_vec_2)
"""
rotation ang
"""
"""
x_mean=float(54.200)
y_mean=float(54.131)
z_increase=4.702
pos_mean=np.array([x_mean,y_mean])
"""
new_base_1=base1.copy()

test_1=base1


layer1=coord.iloc[:1290]
layer2=coord.iloc[1290:]

#arr_1=DataFrame(layer1,columns=['x_coord','y_coord']).to_numpy()
#arr_2=DataFrame(layer2,columns=['x_coord','y_coord']).to_numpy()
#result=[degree(x,pos_mean,y) for x,y in zip(arr_1,arr_2)]
#avg_ang=(sum(result)/len(result))





"""
Test

tmp_4=starter[1935:]
tmp_3=starter[1290:1935]
tmp_2=starter[645:1290]
tmp_1=starter[:645]


test_3=tmp_3.copy()
test_1=tmp_1.copy()

print(test_1['x_coord'])
translate_layer(test_1,trans_vec_1,1)
print(test_1['x_coord'])
arr_1=DataFrame(test_1,columns=['x_coord','y_coord']).to_numpy()
arr_3=DataFrame(test_3,columns=['x_coord','y_coord']).to_numpy()
result=[degree(x,c3_mean[:2],y) for x,y in zip(arr_1,arr_3)]
avg_ang=round(sum(result)/len(result),3)

rotate_layer(test_1,avg_ang, c1_mean, c3_mean, 1)
print(test_1['x_coord'])


End test

"""

#add first layer:
avg_ang=1.458
tmp_4=starter[1935:]
tmp_3=starter[1290:1935]

test_3=tmp_3.copy()
test_4=tmp_4.copy()

names=locals()



translate_layer(test_3,trans_vec_1,1)
translate_layer(test_4,trans_vec_2,1)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,1)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,1)
output_tmp_1=starter.append(test_3)
output_final_1=output_tmp_1.append(test_4)

"""
for i in range(2,3):
    exec('translate_layer(test_3,trans_vec_1,{})'.format(i))
    exec('translate_layer(test_4,trans_vec_2,{})'.format(i))
    exec('rotate_layer(test_3,avg_ang,c1_mean,c3_mean,{})'.format(i))
    exec('rotate_layer(test_4,avg_ang,c2_mean,c4_mean,{})'.format(i))
    exec('output_tmp_{}=output_final_{}.append(test_3)'.format(i,i-1))
    exec('output_final_{}=output_tmp_{}.append(test_3)'.format(i,i))
"""
test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,2)
translate_layer(test_4,trans_vec_2,2)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,2)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,2)
output_tmp_2=output_final_1.append(test_3)
output_final_2=output_tmp_2.append(test_4)


test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,3)
translate_layer(test_4,trans_vec_2,3)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,3)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,3)
output_tmp_3=output_final_2.append(test_3)
output_final_3=output_tmp_3.append(test_4)


test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,4)
translate_layer(test_4,trans_vec_2,4)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,4)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,4)
output_tmp_4=output_final_3.append(test_3)
output_final_4=output_tmp_4.append(test_4)


test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,5)
translate_layer(test_4,trans_vec_2,5)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,5)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,5)
output_tmp_5=output_final_4.append(test_3)
output_final_5=output_tmp_5.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,6)
translate_layer(test_4,trans_vec_2,6)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,6)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,6)
output_tmp_6=output_final_5.append(test_3)
output_final_6=output_tmp_6.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,7)
translate_layer(test_4,trans_vec_2,7)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,7)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,7)
output_tmp_7=output_final_6.append(test_3)
output_final_7=output_tmp_7.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,8)
translate_layer(test_4,trans_vec_2,8)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,8)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,8)
output_tmp_8=output_final_7.append(test_3)
output_final_8=output_tmp_8.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,9)
translate_layer(test_4,trans_vec_2,9)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,9)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,9)
output_tmp_9=output_final_8.append(test_3)
output_final_9=output_tmp_9.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,10)
translate_layer(test_4,trans_vec_2,10)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,10)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,10)
output_tmp_10=output_final_9.append(test_3)
output_final_10=output_tmp_10.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,11)
translate_layer(test_4,trans_vec_2,11)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,11)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,11)
output_tmp_11=output_final_10.append(test_3)
output_final_11=output_tmp_11.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,12)
translate_layer(test_4,trans_vec_2,12)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,12)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,12)
output_tmp_12=output_final_11.append(test_3)
output_final_12=output_tmp_12.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,13)
translate_layer(test_4,trans_vec_2,13)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,13)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,13)
output_tmp_13=output_final_12.append(test_3)
output_final_13=output_tmp_13.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,14)
translate_layer(test_4,trans_vec_2,14)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,14)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,14)
output_tmp_14=output_final_13.append(test_3)
output_final_14=output_tmp_14.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,15)
translate_layer(test_4,trans_vec_2,15)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,15)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,15)
output_tmp_15=output_final_14.append(test_3)
output_final_15=output_tmp_15.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,16)
translate_layer(test_4,trans_vec_2,16)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,16)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,16)
output_tmp_16=output_final_15.append(test_3)
output_final_16=output_tmp_16.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,17)
translate_layer(test_4,trans_vec_2,17)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,17)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,17)
output_tmp_17=output_final_16.append(test_3)
output_final_17=output_tmp_17.append(test_4)

test_3=tmp_3.copy()
test_4=tmp_4.copy()
translate_layer(test_3,trans_vec_1,18)
translate_layer(test_4,trans_vec_2,18)
rotate_layer(test_3,avg_ang,c1_mean,c3_mean,18)
rotate_layer(test_4,avg_ang,c2_mean,c4_mean,18)
output_tmp_18=output_final_17.append(test_3)
output_final_18=output_tmp_18.append(test_4)

"""
for i in range(11,21):
    test_3=tmp_3.copy()
    test_4=tmp_4.copy()
    exec('translate_layer(test_3,trans_vec_1,{})'.format(i))
    exec('translate_layer(test_4,trans_vec_2,{})'.format(i))
    exec('rotate_layer(test_3,avg_ang,c1_mean,c3_mean,{})'.format(i))
    exec('rotate_layer(test_4,avg_ang,c2_mean,c4_mean,{})'.format(i))
    exec('output_tmp_{}=output_final_{}.append(test_3)'.format(i,i-1))
    exec('output_final_{}=output_tmp_{}.append(test_3)'.format(i,i))
"""



write_to_pdb(output_final_18,'elong_20.pdb')







"""
tmp_base=x.copy()
modi_layer(tmp_base, avg_ang, pos_mean, z_increase,1)
output=starter.append(tmp_base)


tmp_base_2=x.copy()
modi_layer(tmp_base_2, avg_ang, pos_mean, z_increase,2)
output_2=output.append(tmp_base_2)



tmp_base_3=x.copy()
modi_layer(tmp_base_3, avg_ang, pos_mean, z_increase,3)
output_3=output_2.append(tmp_base_3)

tmp_base_4=x.copy()
modi_layer(tmp_base_4, avg_ang, pos_mean, z_increase,4)
output_4=output_3.append(tmp_base_4)

tmp_base_5=x.copy()
modi_layer(tmp_base_5, avg_ang, pos_mean, z_increase,5)
output_5=output_4.append(tmp_base_5)

tmp_base_6=x.copy()
modi_layer(tmp_base_6, avg_ang, pos_mean, z_increase,6)
output_6=output_5.append(tmp_base_6)

tmp_base_7=x.copy()
modi_layer(tmp_base_7, avg_ang, pos_mean, z_increase,7)
output_7=output_6.append(tmp_base_7)

tmp_base_8=x.copy()
modi_layer(tmp_base_8, avg_ang, pos_mean, z_increase,8)
output_8=output_7.append(tmp_base_8)

tmp_base_9=x.copy()
modi_layer(tmp_base_9, avg_ang, pos_mean, z_increase,9)
output_9=output_8.append(tmp_base_9)

tmp_base_10=x.copy()
modi_layer(tmp_base_10, avg_ang, pos_mean, z_increase,10)
output_10=output_9.append(tmp_base_10)

tmp_base_11=x.copy()
modi_layer(tmp_base_11, avg_ang, pos_mean, z_increase,11)
output_11=output_10.append(tmp_base_11)

tmp_base_12=x.copy()
modi_layer(tmp_base_12, avg_ang, pos_mean, z_increase,12)
output_12=output_11.append(tmp_base_12)

tmp_base_13=x.copy()
modi_layer(tmp_base_13, avg_ang, pos_mean, z_increase,13)
output_13=output_12.append(tmp_base_13)

tmp_base_14=x.copy()
modi_layer(tmp_base_14, avg_ang, pos_mean, z_increase,14)
output_14=output_13.append(tmp_base_14)

tmp_base_15=x.copy()
modi_layer(tmp_base_15, avg_ang, pos_mean, z_increase,15)
output_15=output_14.append(tmp_base_15)

tmp_base_16=x.copy()
modi_layer(tmp_base_16, avg_ang, pos_mean, z_increase,16)
output_16=output_15.append(tmp_base_16)

tmp_base_17=x.copy()
modi_layer(tmp_base_17, avg_ang, pos_mean, z_increase,17)
output_17=output_16.append(tmp_base_17)

tmp_base_18=x.copy()
modi_layer(tmp_base_18, avg_ang, pos_mean, z_increase,18)
output_18=output_17.append(tmp_base_18)

tmp_base_19=x.copy()
modi_layer(tmp_base_19, avg_ang, pos_mean, z_increase,19)
output_19=output_18.append(tmp_base_19)

tmp_base_20=x.copy()
modi_layer(tmp_base_19, avg_ang, pos_mean, z_increase,20)
output_20=output_19.append(tmp_base_20)

write_to_pdb(output_20,'elong_20.pdb')
"""


'''
modi_layer(new_base_1,avg_ang,pos_mean,z_increase)

output=starter.append(new_base_1,ignore_index=True)

write_to_pdb(output,'3layers.pdb')'''