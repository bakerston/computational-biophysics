if [ $# -eq 0 ]
then   
	echo "usage: ./find_charge.sh pdb"
	exit
else
	new_end=$1

fi
echo $1
CHAIN_ORI=$( grep -o "TER" $1 | wc -l )
echo $CHAIN_ORI
rm -rf charge.list

ARG=$( awk 'BEGIN {n=0} $3=="CA" && $4=="ARG" {n++} END{print n}' $1 )
LYS=$( awk 'BEGIN {n=0} $3=="CA" && $4=="LYS" {n++} END{print n}' $1 ) 
ASP=$( awk 'BEGIN {n=0} $3=="CA" && $4=="ASP" {n++} END{print n}' $1 )  
GLU=$( awk 'BEGIN {n=0} $3=="CA" && $4=="GLU" {n++} END{print n}' $1 )

NA=$( awk 'BEGIN {n=0} $3=="NA" {n++} END{print n}' $1 ) 
CL=$( awk 'BEGIN {n=0} $3=="CL" {n++} END{print n}' $1 ) 

ANS=$[$ARG+$LYS+$NA-$ASP-$GLU-$CL]
echo $ANS

if [ $ANS -lt 0 ]
then
	echo "NEG"
	NEG= $[$ASP+$GLU+$CL-$ARG-$LYS-$NA]	
	for i in $(seq 1 $NEG)
	do	
		cat NA.pdb >> $1
       		echo "$[i + $CHAIN_ORI].1.NA   1" >> charge.list
	done
elif [ $ANS -gt 0 ]
then
	echo "POS"
	for i in $(seq 1 $ANS)
	do
		cat CL.pdb >> $1
		echo "$[i + $CHAIN_ORI].1.CL   -1" >> charge.list
	done
else
	echo "NEU"
fi

	#old_end=`head -n 2 run.pbs | tail -n 1 | awk -F '[-]' '{print $NF}'`

#task_num=`head -n 7 run.pbs | tail -n 1 | awk -F '[=:]' '{print $2}'`


#if [ $old_end -ge $new_end ]
#then
#	echo "retype new_end"
#	exit
#else
#	cur_step=$[ $new_end - $old_end ]
#fi
#echo $[cur_step]

#sed  -i "2s/$[old_end]/$[new_end]/" run.pbs
#sed  -i "s/-$[old_end]/-$[new_end]/g" equi.task

#old_step=`tail -n 1 equi.task | awk -F '[ ]' '{print $NF}'`


#sed  -i "s/$[old_step]/$[$cur_step * 1000000]/g" equi.task

#rm -rf tmp.jobs

#for i in $(seq 1 $[$task_num - 1]);do echo "cd" `pwd`"/$i;~/dmd-1.1/apps/xDMD.linux -p Para_out -s equi-$[$old_end]m.dmd_restart -c Constr_out -i ../equi.task";done > tmp.jobs

