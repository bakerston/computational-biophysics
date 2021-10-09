# coding=utf-8
import math
import numpy as np
import sys
from scipy.spatial.distance import cdist
"""
读取所有CA坐标
81 chains * 37 residues * 100 snaps
CA.1.1 CA.2.1 为fibril方向矢量

CA.81.* 算平均坐标为monomer位置
"""

#dist = numpy.linalg.norm(a-b)
#convert input string to list.
def str2list(input_str):
   return input_str.split()

def t(p, q, r):
    x = p-q
    return np.dot(r-q, x)/np.dot(x, x)

#get distance from point r to line section given by p, q
def dis(p, q, r):
    return np.linalg.norm(t(p, q, r)*(p-q)+q-r)

def get_com(data_in, data_len):
    ans=[]

    x_list=[x[0] for x in data_in]
    y_list=[x[1] for x in data_in]
    z_list=[x[2] for x in data_in]

    ans.append(sum(x_list)/data_len)
    ans.append(sum(y_list)/data_len)
    ans.append(sum(z_list)/data_len)

    return ans

def mid_p(p_1,p_2):
    return list(map(lambda x,y:(x+y)/2,p_1,p_2))

def dis_along_f(p,q,m,x):
    perp=dis(p,q,x)
    full=np.linalg.norm(x-m)
    return math.sqrt(full*full-perp*perp)

#setup

f=open("3-3m.pdb")
line=f.readline()
tmp=[]
res_num,ori_mark=37,1


#read fibril vector for 1st time
while ori_mark<=1000:
    if line[0] in ["E","T"]:
        pass
    else:
        if line[13:15] != "CA":
            pass
        else:
            cur=str2list(line)
            tmp.append(list(map(float,line[31:55].split())))
            ori_mark+=1
    line=f.readline()
f.close()

m_18=get_com(tmp[629:666],res_num)
m_19=get_com(tmp[666:703],res_num)
m_21=get_com(tmp[740:777],res_num)
m_22=get_com(tmp[777:814],res_num)

end_a=np.array(mid_p(m_18,m_19))
end_b=np.array(mid_p(m_21,m_22))

center=np.array(mid_p(end_a,end_b))

print(end_b,end_a,center)
tmp,mono_pos=[],[]

f=open("3-3m.pdb")
line_2=f.readline()
sec_mark=0


while line_2:
    if line_2[0] in ["E","T"]:
        pass

    else:
        if line_2[13:15] != "CA":
            pass

        else:
            if sec_mark < 2960:
                sec_mark += 1

            elif sec_mark < 2997:
                tmp.append(list(map(float,line_2[31:55].split())))
                sec_mark += 1

            else:
                sec_mark = 1
                mono = get_com(tmp, 37)
                tmp = []
                mono_pos.append(mono)
    line_2 = f.readline()
f.close()

mono_pos.append(get_com(tmp,37))

mono=[np.array(x) for x  in mono_pos]

ans=[dis_along_f(end_a,end_b,center,x) for x in mono]




f=open("long_ca.dat","w+")
for i in range(len(ans)):
    f.write(str(ans[i]))
    f.write("\n")
f.close()