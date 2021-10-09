#helix beta turn bend coil
for i in $(seq 1 50)
	do cd ${i}/dssp/
	rm -f all_ss.dat
	for j in {2..10..2}
		do
		sed -n '12,211p' ss_${j}.dat>> all_ss.dat
		
        done
	awk '{print $2}' all_ss.dat | cut -d'!' -f21 >final.dat
	
	echo ${i}
	awk -F'|' '{a=gsub(/H/,"")+gsub(/G/,"")+gsub(/I/,"");b=gsub(/B/,"")+gsub(/E/,"");c=gsub(/T/,"");d=gsub(/S/,"")+gsub(/C/,"")}{printf("%12s\t",a/37)}{printf("%12s\t",b/37)}{printf("%12s\t", c/37)} {printf("%12s\t", d/37)}' final.dat > freq.dat
	cd ../..
done	
