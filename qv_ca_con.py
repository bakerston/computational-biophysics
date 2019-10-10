import math

#convert input string to list.
def str2list(input_str):
   return input_str.split()

#calculate 3D distance of given point set.
def dis(x,y):
   return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)

#exist contacts among two particular residues
def is_cont(res_1,res_2,dis_cut):
    return 1 if dis(res_1.ca,res_2.ca) <= dis_cut else 0

#x-mo,y-fi
def get_q_value(data_in,dis_cut,res_num):
    tmp = [[0 for m in range(res_num)] for n in range(res_num)]
    fi_list=data_in[:740]
    mo_list=data_in[740:]
    for fi in range(len(fi_list)):
        for mo in range(len(mo_list)):
            res_1=residue(fi_list[fi][0],fi_list[fi][1])
            res_2=residue(mo_list[mo][0]%37,mo_list[mo][1])
            tmp[(fi%37)-1][mo]+=is_cont(res_1,res_2,dis_cut)
    return tmp

class residue:
    def __init__(self, res_index, res_ca):
        self.index = res_index
        self.ca=res_ca

def cal_q_value(matrix_in,eff_range):
    ans=0
    tmp=0
    for i in range(len(matrix_in)):
        for j in range(i-eff_range,i+eff_range+1):
            try:
                tmp+=matrix_in[i][j]
            except IndexError:
                continue
        if tmp>0:
            ans+=1
        tmp=0
    return ans/len(matrix_in)


#setup
dis_cut=8.5
res_num=37
eff_range=2

snap_mark=0
chain_mark=0

cur_snap=[]

data_in=[]
cur_data=[]
q_list=[]
i=0

f=open("10-200steps.pdb")
line=f.readline()

while line:
    cur_list=str2list(line)
    if cur_list[0]=="TER":
        line=f.readline()
        continue
    if cur_list[0]=="ENDMDL":
        snap_mark+=1
        print(snap_mark)

        tmp=get_q_value(data_in,dis_cut,res_num)

        out_q=cal_q_value(tmp,eff_range)

        q_list.append(out_q)

        data_in=[]

    else:
        if  cur_list[2]=="CA":
            cur_data.append(int(cur_list[5]))
            cur_data.append(list(map(float, cur_list[6:9])))
            data_in.append(cur_data)
            cur_data = []
    line = f.readline()

#print(ans)
f.close()
res=0



f=open("qv_ca_snap.dat","w+")
for i in range(len(q_list)):
    f.write(str(q_list[i]))
    f.write("\n")
f.write("max:")
f.write("\t")
f.write(str(max(q_list)))

f.close()