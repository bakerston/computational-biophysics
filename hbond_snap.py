import math
from collections import Counter
import sys

#convert input string to list.
def str2list(input_str):
   return input_str.split()

#calculate 3D angle of target point set.
def angle(x,y,z):
   vec_1=list(map(lambda i,j: i-j, x,y))
   vec_2=list(map(lambda i,j: i-j, z,y))
   ans=sum(map(lambda x,y: x*y,vec_1,vec_2))/(math.sqrt(sum(map(lambda x:x*x,vec_1)))*math.sqrt(sum(map(lambda y:y*y,vec_2))))
   return abs(math.acos(ans)*180/math.pi)

#calculate 3D distance of given points.
def distance(x,y):
   return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)

#output all legit hbond pairs
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
    return ans

def get_num_pairs(mono_list,chain_list,cut_off,angle_limit):
    ans=0
    for mono in range(len(mono_list)):
        for chain in range(len(chain_list)):
            if distance(mono_list[mono][0],chain_list[chain][1])<=cut_off and angle(mono_list[mono][0],mono_list[mono][2],chain_list[chain][1]) >= angle_limit:
                ans+=1
    #mono as acceptor
    for mono in range(len(mono_list)):
        for chain in range(len(chain_list)):
            if distance(mono_list[mono][1], chain_list[chain][0]) <= cut_off and angle(chain_list[chain][0],chain_list[chain][2],mono_list[mono][1]) >= angle_limit:
                ans+=1
    return ans


#count snap number in the given file.
def file_snaps(file,key_word):
    line=file.readline()
    mark=0
    while line:
        if str2list(line)[0]==key_word:
            mark+=1
        line=file.readline()
    return mark

#calculate and output frequency number of hbond on each residue.
def sum_freq(data_in,save_name):
    mono = []
    tmp_chain = []
    for x in range(len(data_in)):
        for y in range(len(data_in[x])):
            mono.append(data_in[x][y][0])
            tmp_chain.append(data_in[x][y][1])
    chain = [x % 37 for x in tmp_chain]
    mono_count = Counter(mono)
    mono_comb = list(zip(mono_count.keys(), mono_count.values()))
    print(mono_comb)
    f=open(save_name,"w+")
    ans=list(map(lambda x: list(x),mono_comb))
    for i in range(len(ans)):
        for j in range(len(ans[i])):
            f.write(str(ans[i][j]))
            f.write('\t')
        f.write('\n')
    f.close()

#output
def sored_count_freq(data_in):
    mono=[]
    tmp_chain=[]
    for x in range(len(data_in)):
        for y in range(len(data_in[x])):
            mono.append(data_in[x][y][0])
            tmp_chain.append(data_in[x][y][1])
    chain=[x%37 for x in tmp_chain]
    mono_count=Counter(mono)
    chain_count=Counter(chain)

    mono_comb = list(zip(mono_count.keys(), mono_count.values()))
    mono_comb.sort(key=lambda x: x[0])

    chain_comb= list(zip(chain_count.keys(),chain_count.values()))
    chain_comb.sort(key=lambda x: x[0])
    return [mono_comb, chain_comb]

#output all pairs
#def print_pairs(data_in):
#    f=open("pair.dat","w+")
#    for i in data_in:
#        np.savetxt("pair.dat",data_in[i])
#    f.close()
#    for i in data_in:
#        f.write("%s\n\t" %i)
#    f.close()


#output hbond numbers for each snap.
def print_num_pairs(data_in,save_name):
    f=open(save_name,"w+")
    for i in data_in:
        alen=len(i)
        f.write("%d\n" %alen)
    f.close()

def print_num_pairs_lite(data_in,save_name):
    f = open(save_name, "w+")
    for i in data_in:
        f.write("%d\n" % i)
    f.close()

def main_read(load_name, all_chain_snap,all_mono_snap):
    f=open(load_name)
    snap_count=file_snaps(f,"ENDMDL")
    print(snap_count)
#    cut_off=3.5
#    ang_limit=120

#    all_chain_snap=[]
#    all_mono_snap=[]

    cur_chain=[]
    cur_mono=[]
    single_line=[]

    res_mark=0
    snap_mark=0

    while snap_mark<snap_count:
        cur_list=str2list(f.readline())
        if not cur_list:
            continue
        elif cur_list[0]=="TER" or cur_list[0]=="ENDMDL":
            if cur_list[0]=="ENDMDL":
                all_chain_snap.append(cur_chain)
                cur_chain=[]
                all_mono_snap.append(cur_mono)
                cur_mono=[]
            else:
                res_mark+=1
                print(res_mark)
                if res_mark>20:
                    res_mark=0
                    snap_mark+=1

        else:
            if cur_list[2]=="N":
                single_line.append(list(map(float, cur_list[6:9])))
            if cur_list[2]=="O":
                single_line.append(list(map(float, cur_list[6:9])))
            if cur_list[2]=="HN":
                single_line.append(list(map(float, cur_list[6:9])))
                if res_mark<20:
                    cur_chain.append(single_line)
                    single_line=[]
                else:
                    cur_mono.append(single_line)
                    single_line=[]
        all_chain_snap.append(cur_chain)
        cur_chain=[]
        all_mono_snap.append(cur_mono)
        cur_mono=[]
        print(all_mono_snap)
    return all_chain_snap,all_mono_snap,snap_count



def get_all_pairs(all_chain_snap,all_mono_snap,cut_off,ang_limit,snap_count):
    ans=[]
    for snap in range(snap_count):
        tmp=get_pairs(all_mono_snap[snap],all_chain_snap[snap],cut_off,ang_limit)
        ans.append(tmp)
    return ans

def get_num_all_pairs(all_chain_snap, all_mono_snap, cut_off, ang_limit, snap_count):
    ans_num=[]
    for snap in range(snap_count):
        tmp=get_num_pairs(all_mono_snap[snap],all_chain_snap[snap],cut_off,ang_limit)
        ans_num.append(tmp)
    return ans_num

main_read("1-4m.pdb",[],[])


"""

usage=sys.argv[:]

if len(usage)==4:
   # if type(usage[1])==type(usage[5])==int and type(usage[2])==type(usage[3])==str and type(usage[4])==float:
    if usage[1]==1:
        #chain,mono,snap_count=main_read(usage[2])
        #ans=get_all_pairs(chain,mono,3.5,120,snap_count)
        pass

    elif int(usage[1])==2:
        print("in progress 2")
        chain, mono, snap_count = main_read(usage[2])
        ans = get_all_pairs(chain, mono, 3.5,120,snap_count)
        sum_freq(ans,usage[3])

    elif int(usage[1])==3:
        print("in progress 3")
        chain, mono, snap_count = main_read(usage[2])
        print("read")
        ans = get_num_all_pairs(chain, mono, 3.5,120,snap_count)
        print("ans")
        print_num_pairs_lite(ans,usage[3])
        print("done")

    else:
        print("usage: python hbond_snap.py job_type l_n s_n ")
else:
    print("usage: python hbond_snap.py job_type l_n s_n ")
"""
