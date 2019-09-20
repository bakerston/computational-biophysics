import math
from collections import Counter
import sys

#convert input string to list.
def str2list(input_str):
   return input_str.split()

#calculate 3D angle of target point set.
def ang(x,y,z):
   vec_1=list(map(lambda i,j: i-j, x,y))
   vec_2=list(map(lambda i,j: i-j, z,y))
   ans=sum(map(lambda x,y: x*y,vec_1,vec_2))/(math.sqrt(sum(map(lambda x:x*x,vec_1)))*math.sqrt(sum(map(lambda y:y*y,vec_2))))
   return abs(math.acos(ans)*180/math.pi)

#calculate 3D distance of given point set.
def dis(x,y):
   return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)


"""#record all legit hbond pairs in one snap
def get_pairs(mono_list,chain_list,cut_off,angle_limit):
    ans=[]
    #mono as donor
    for mono in range(len(mono_list)):
        for chain in range(len(chain_list)):
            if distance(mono_list[mono][0],chain_list[chain][1])<=cut_off and angle(mono_list[mono][0],mono_list[mono][2],chain_list[chain][1]) >= angle_limit:
                ans.append([mono,chain])
    #mono as acceptor
    for mono in range(len(mono_list)):
        for chain in range(len(chain_list)):
            if distance(mono_list[mono][1], chain_list[chain][0]) <= cut_off and angle(chain_list[chain][0],chain_list[chain][2],mono_list[mono][1]) >= angle_limit:
                ans.append([mono,chain])
    return ans"""

#number of hbonds among two particular residues
def num_hbond(res_1,res_2,dis_cut,ang_cut):
    ans=0
    if dis(res_1.n,res_2.o)<=dis_cut and ang(res_1.n,res_1.h,res_2.o)>=ang_cut:
        ans+=1
    if dis(res_2.n,res_1.o)<=dis_cut and ang(res_2.n,res_2.h,res_1.o)>=ang_cut:
        ans+=1
    return ans

#x-mo,y-fi
def get_q_value(data_in,dis_cut,ang_cut,res_num,eff_range):
    tmp = [[0 for m in range(res_num)] for n in range(res_num)]
    fi_list=data_in[:740]
    mo_list=data_in[740:]
    for fi in range(len(fi_list)):
        for mo in range(len(mo_list)):
            res_1=residue(fi_list[fi][0],fi_list[fi][1],fi_list[fi][2],fi_list[fi][3])
            res_2=residue(mo_list[mo][0]%37,mo_list[mo][1],mo_list[mo][2],mo_list[mo][3])

            tmp[fi%37-1][mo]+=num_hbond(res_1,res_2,dis_cut,ang_cut)



    return tmp

#def count_hbond(res_1,res_2,dis_cut,ang_cut):
#    if is_hbond(res_1,res_2,dis_cut,ang_cut):
#        return

class residue:
    def __init__(self, res_index, n_pos,o_pos,h_pos):
        self.index = res_index
        self.n = n_pos
        self.h = h_pos
        self.o = o_pos

def file_snaps(file,key_word):
    line=file.readline()
    mark=0
    while line:
        if str2list(line)[0]==key_word:
            mark+=1
        line=file.readline()
    return mark


#setup
dis_cut=3.5
ang_cut=120
res_num=37
eff_range=0
ans=[[0 for x in range(res_num)] for y in range(res_num)]



#snap_count=file_snaps(f,"ENDMDL")
snap_mark=0
chain_mark=0

cur_snap=[]

#for each snap within "ENDMDL"
#读取数据，储存至两个数组，fi和mo


#for fi_chain i:
#   for fi_res j:
#       for mo_res k:
#           if is_hb(j,k):
#               tmp_list[j][k]+=num_hb(j,k)
#
#if cur_list="ENDMDL"
#   清空fi和mo
#   将tmp_list写入ans

data_in=[]
cur_data=[]

i=0
print(ans)
f=open("100steps.pdb")
line=f.readline()
while line:
    cur_list=str2list(line)
    if cur_list[0]=="TER":
        line=f.readline()
        continue
    if cur_list[0]=="ENDMDL":
        snap_mark+=1
        print(snap_mark)
        tmp=get_q_value(data_in,dis_cut,ang_cut,res_num,eff_range)
        for i in range(len(ans)):
            for j in range(len(ans[0])):
                ans[i][j]+=tmp[i][j]
        data_in=[]

    else:
        if cur_list[2] == "N":
            cur_data.append(int(cur_list[5]))
            cur_data.append(list(map(float, cur_list[6:9])))
        if cur_list[2] == "O":
            cur_data.append(list(map(float, cur_list[6:9])))
        if cur_list[2] == "HN":
            cur_data.append(list(map(float, cur_list[6:9])))
            data_in.append(cur_data)
            cur_data = []
    line = f.readline()
print(ans)
#print(ans)
f.close()
res=0
for i in range(len(ans)):
    for j in range(len(ans[0])):
        res+=ans[i][j]
print(res)

f=open("xyz-q.dat","w+")
for i in range(len(ans)):
    for j in range(len(ans[0])):
        f.write(str(i))
        f.write("\t")
        f.write(str(j))
        f.write("\t")
        f.write(str(ans[i][j]))
        f.write("\n")
f.close()