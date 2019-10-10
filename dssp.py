# coding=utf-8

def str2list(input_str):
   return input_str.split()

alist=["H","G","I"]
blist=["B","E"]
tlist=["T"]
clist=["S","C"]

snap_num=200
file_num=50
resi_num=37

a_outlist=[[0 for n in range(resi_num)] for i in range(file_num)]
b_outlist=[[0 for n in range(resi_num)] for j in range(file_num)]
prelist=[0 for n in range(resi_num)]
t_outlist=[[0 for n in range(resi_num)] for i in range(file_num)]
c_outlist=[[0 for n in range(resi_num)] for j in range(file_num)]


for i in range(1,file_num+1):

    f=open("10-"+str(i)+".dat")
    a_fo=open("ahelix-"+str(i)+".dat","w+")
    b_fo=open("bsheet-"+str(i)+".dat","w+")
    t_fo=open("turn-"+str(i)+".dat","w+")
    c_fo=open("coil-"+str(i)+".dat","w+")
    line=f.readline()
    while line:
        tmp=str2list(line)

        if tmp[1] in alist:
            a_outlist[i-1][int(tmp[0])-1]+=int(tmp[2])
        if tmp[1] in blist:
            b_outlist[i-1][int(tmp[0])-1]+=int(tmp[2])
        if tmp[1] in tlist:
            t_outlist[i-1][int(tmp[0])-1]+=int(tmp[2])
        if tmp[1] in clist:
            c_outlist[i-1][int(tmp[0])-1]+=int(tmp[2])
        line=f.readline()
    print(i)
    for out in range(len(a_outlist[i-1])):
        a_fo.write(str(a_outlist[i-1][out]/snap_num))
        a_fo.write("\t")
    a_fo.write("\n")
    a_fo.close()

    for out in range(len(b_outlist[i-1])):
        b_fo.write(str(b_outlist[i-1][out]/snap_num))
        b_fo.write("\t")
    b_fo.write("\n")
    b_fo.close()

    for out in range(len(t_outlist[i-1])):
        t_fo.write(str(t_outlist[i-1][out]/snap_num))
        t_fo.write("\t")
    t_fo.write("\n")
    t_fo.close()

    for out in range(len(c_outlist[i-1])):
        c_fo.write(str(c_outlist[i-1][out]/snap_num))
        c_fo.write("\t")
    c_fo.write("\n")
    c_fo.close()
    f.close()
