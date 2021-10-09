if [ $# -eq 0 ]
then
	echo "usage: ./shell func"
	exit
#elif	[ $1 = "hbond" ]
#	start_time=` sed -n '2p' try.data | awk '{print $1;}'`                              end_time=` tail -n 1 try.data | awk '{print $1;}'`                                                                                                                                                                                                          end=`echo "($end_time-$start_time)/100 + 1" | bc` 
#	n=` tail -n 1 try.data | awk '{print $1;}'`                                         echo $n                                                                                                                                                                                                                                                     end=`echo "$n/100" | bc`                                                            echo $end        
#then
#	for i in $(seq 1 50);
#	do 
#		cd  $i
#	start_time=` sed -n '2p'  | awk '{print $1;}'`
#	end_time=` tail -n 1 try.data | awk '{print $1;}'` 
#	echo "cd" `pwd`"/$i; /zfs/biophys/czhang6/medusa-1.1/src/apps/complex_M2P.linux -P ~/medusa-1.1/parameter/ -I hif-abd-$i.pdb -M equi-2m.dmd_movie -o ${i}-200steps.pdb -F 1:100:100 ";done > hbond.jobs 
elif    [ $1 = "hbond" ]
then   
   
	mkdir hbond	
	for i in 2 4 6 8 10
	do 
		rm -rf hbond-${i}m.task 
	   	echo "HBOND_BB   hb-abinter-${i}m    1.*.*    2.*.*" > hbond-${i}m.task
	   	echo "HBOND_BB   hb-intro-1-${i}m    1.*.*    3-22.*.*" >> hbond-${i}m.task 
	   	echo "HBOND_BB   hb-intro-2-${i}m    2.*.*    3-22.*.*" >> hbond-${i}m.task
           	echo "HBOND_BB   hb-intro-d-${i}m    1-2.*.*  3-22.*.*" >> hbond-${i}m.task
	   	cp hbond-${i}m.task hbond/
		echo task-$i
		cp run.pbs hbond-${i}m.pbs
		sed -i "2c #PBS -N hbond-${i}m" hbond-${i}m.pbs
		sed -i "s/dmd-10/hbond-${i}m/" hbond-${i}m.pbs 
		echo pbs-$i
   	done

	for j in 2 4 6 8 10
	do
	for i in $(seq 1 50)
	do
		echo "cd" `pwd`"/$i/hbond; ~/medusa-1.1/src/dev-tools/complex_analysis.linux -P ~/medusa-1.1/parameter/ -I ../hif-abd-$i.pdb -M ../equi-${j}m.dmd_movie -t hbond-${j}m.task -T /zfs/biophys/czhang6/newstart/mol2/NACL.topparam -F 1:20000:1 ";done > hbond-${j}m.jobs 
	echo jobs-$j
	done
	for i in $(seq 1 50)
		do 
		cd $i
		rm -rf hbond	
		cp -r ../hbond .
		cd ..
		echo cp-$i
	done

	elif	[ $1 = "contact" ]
then       
	mkdir contact	
	for i in 2 4 6 8 10
	do 
		rm -rf contact-${i}m.task 
	   	echo "CONTACTS   cont-abinter-${i}m    1.*.*    2.*.*"      6.5 >  contact-${i}m.task
	   	echo "CONTACTS   cont-intro-1-${i}m    1.*.*    3-22.*.*"   6.5 >> contact-${i}m.task 
	   	echo "CONTACTS   cont-intro-2-${i}m    2.*.*    3-22.*.*"   6.5 >> contact-${i}m.task
           	echo "CONTACTS   cont-intro-d-${i}m    1-2.*.*  3-22.*.*"   6.5 >> contact-${i}m.task
	   	cp contact-${i}m.task contact/
		echo task-$i
		cp run.pbs contact-${i}m.pbs
		sed -i "2c #PBS -N contact-${i}m" contact-${i}m.pbs
		sed -i "s/dmd-12/contact-${i}m/" contact-${i}m.pbs 
		echo pbs-$i
   	done

	for j in 2 4 6 8 10
	do
	for i in $(seq 1 50);do echo "cd" `pwd`"/$i/contact; ~/medusa-1.1/src/dev-tools/complex_analysis.linux  -P ~/medusa-1.1/parameter/ -I ../hif-abd-$i.pdb -M ../equi-${j}m.dmd_movie -t contact-${j}m.task -T /zfs/biophys/czhang6/newstart/mol2/NACL.topparam -F 1:20000:1 ";done > contact-${j}m.jobs 
	done

	for i in $(seq 1 50)
		do 
		cd $i
	
		cp -r ../contact .
		cd ..
	done

	elif	[ $1 = "dssp" ]
then
	for i in $(seq 1 50);do echo "cd" `pwd`"/$i; ~/medusa-1.1/src/apps/complex_M2P.linux -P ~/medusa-1.1/parameter/ -I hif-abd-$i.pdb -M equi-1m.dmd_movie -o ${i}-100steps.pdb -F 1:100:100 ";done > dssp.jobs 

	
	elif	[ $1 = "m2p" ]

then
	for j in 2 4 6 8 10 
	do
	for i in $(seq 1 20);do echo "cd" `pwd`"/$i; ~/medusa-1.1/src/apps/complex_M2P.linux -P ~/medusa-1.1/parameter/ -I abf-abd-$i.pdb -M equi-${j}m.dmd_movie -o ${i}-${j}00steps.pdb -T /zfs/biophys/czhang6/newstart/mol2/NACL.topparam -F 1:200:100 ";done > m2p-${j}m.jobs 
	cp run.pbs m2p-${j}m.pbs
	sed -i "2c #PBS -N m2p-${j}m" m2p-${j}m.pbs
	sed -i "s/dmd/m2p-${j}m/" m2p-${j}m.pbs 
	echo pbs-$i
done

	elif	[ $1 = "dmd" ]
then
	for i in $(seq 1 20);do echo "cd" `pwd`"/$i;~/dmd-1.1/apps/xDMD.linux -p Para_out -s equi-10m.dmd_restart -c Constr_out -i ../equi.task";done > dmd.jobs 
else
	echo "input valid func"
fi
