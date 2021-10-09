for j in $(seq 1 20)
    do
	mkdir $j
	/zfs/biophys/czhang6/medusa-1.1/src/apps/randomize_complexSelect.linux    /zfs/biophys/czhang6/medusa-1.1/parameter/  abf-abd.pdb 150  $[j*123] ${j}/abf-abd-${j}.pdb 15	3-22 /zfs/biophys/czhang6/newstart/mol2/NACL.topparam 
    done
 

    


for j in $(seq 1 1)
    do
	cd $j
	~/medusa-1.1/src/dev-tools/genTriplet.linux   ~/medusa-1.1/parameter/   abf-abd-${j}.pdb > pdbConstr
	~/medusa-1.1/src/apps/genESC.linux -P   ~/medusa-1.1/parameter/  -I abf-abd-${j}.pdb -T /zfs/biophys/czhang6/newstart/mol2/NACL.topparam -L /zfs/biophys/czhang6/newstart/mol2/charge.list > escConstr
	~/medusa-1.1/src/apps/findSS.linux  abf-abd-${j}.pdb > disConstr
	cat disConstr pdbConstr escConstr > Constr_in
	cd ..
done
	echo "Static 3-22.*.*"  > tmp
        cat tmp 1/Constr_in > Constr_in 

for j in $(seq 1  1)
	do
	cp Constr_in $j/
	cd $j
	~/medusa-1.1/src/apps/complex.linux -P   ~/medusa-1.1/parameter/ -I abf-abd-${j}.pdb -D 150 -p Para_out -s Stat_out -T /zfs/biophys/czhang6/newstart/mol2/NACL.topparam -C Constr_in -c Constr_out
	cd ..
    done
