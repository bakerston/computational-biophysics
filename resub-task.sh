if [ $# -eq 0 ]
then
	echo "usage: ./shell new_end(M)"
	exit
else
	new_end=$1

fi

old_end=`head -n 2 run.pbs | tail -n 1 | awk -F '[-]' '{print $NF}'`

task_num=`head -n 7 run.pbs | tail -n 1 | awk -F '[=:]' '{print $2}'`


if [ $old_end -ge $new_end ]
then
	echo "retype new_end"
	exit
else
	cur_step=$[ $new_end - $old_end ]
fi
echo $[cur_step]

sed  -i "2s/$[old_end]/$[new_end]/" run.pbs
sed  -i "s/-$[old_end]/-$[new_end]/g" equi.task

old_step=`tail -n 1 equi.task | awk -F '[ ]' '{print $NF}'`


sed  -i "s/$[old_step]/$[$cur_step * 1000000]/g" equi.task

rm -rf tmp.jobs

for i in $(seq 1 $[$task_num - 1]);do echo "cd" `pwd`"/$i;~/dmd-1.1/apps/xDMD.linux -p Para_out -s equi-$[$old_end]m.dmd_restart -c Constr_out -i ../equi.task";done > tmp.jobs

