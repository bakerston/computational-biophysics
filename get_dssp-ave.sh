#helix beta turn bend
for i in ahelix bsheet coil turn
	do cd ${i}
        rm -rf ${i}-ave.dat	
	for j in {8,10}
		do cd ${j}
		cat ${i}-*.dat > ${i}-total-${j}
        	awk '{for(i=1; i<=NF; i++) {a[i]+=$i; if($i!="") b[i]++}}; END {for(i=1; i<=NF; i++) printf "%s%s", a[i]/b[i], (i==NF?ORS:OFS)}' ${i}-total-${j} > tmp  
		#awk '{for(i=1;i<=NF;++i){a[i]+=$i}}1;END{print "-avg-";for(x in a){printf "%s",a[x]/NR" "}print ""}' ${i}-total-${j} > tmp
		tail -n 1 tmp> ${i}-${j}-ave.dat
		echo $i-$j
		cd ..
	done
	cat 8/${i}-8-ave.dat 10/${i}-10-ave.dat> ${i}-ave.dat
	cd ..
done	
