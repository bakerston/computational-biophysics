	#helix beta turn bend coil
for i in $(seq 1 50)
	do
	rm -f all_ss.dat
	rm -f final.dat
	rm -f 6m8m_ss.dat
	rm -f 8m10m_ss.dat
	rm -f 8-${i}.dat	
	rm -f 10-${i}.dat
	cp ../../${i}/dssp/ss_8.dat .
	sed -n '12,211p' ss_8.dat>> 6m8m_ss.dat
	wc -l 6m8m_ss.dat	
  
	awk '{print $2}' 6m8m_ss.dat | cut -d'!' -f21 > f21.dat
	
	cat f21.dat | awk '{split($1, chars, "") ;for (i=1;i<=length($1);i++){printf("%s\t",chars[i])}print("\n")} '   > split-f21.dat 
	awk 'NF>0' split-f21.dat > final-f21.dat


	for res in $(seq 1 37)
		 do awk -v I=$res -F '\t' '{sum[$I]++}END{for (i in sum) print I "\t" i "\t" sum[i]}' final-f21.dat >> 8-${i}.dat
	done
	echo 8-${i}
	

	cp ../../${i}/dssp/ss_10.dat .
	sed -n '12,211p' ss_10.dat>> 8m10m_ss.dat
	wc -l 8m10m_ss.dat	
  
	awk '{print $2}' 8m10m_ss.dat | cut -d'!' -f21 > f21.dat
	
	cat f21.dat | awk '{split($1, chars, "") ;for (i=1;i<=length($1);i++){printf("%s\t",chars[i])}print("\n")} '   > split-f21.dat 
	awk 'NF>0' split-f21.dat > final-f21.dat


	for res in $(seq 1 37)
		 do awk -v I=$res -F '\t' '{sum[$I]++}END{for (i in sum) print I "\t" i "\t" sum[i]}' final-f21.dat >> 10-${i}.dat
	done
	echo 10-${i}
done


	#awk -F'|' '{a=gsub(/H/,"")+gsub(/G/,"")+gsub(/I/,"");b=gsub(/B/,"")+gsub(/E/,"");c=gsub(/T/,"");d=gsub(/S/,"");e=gsub(/C/,"")}{printf("%12s\t",a/37)}{printf("%12s\t",b/37)}{printf("%12s\t", c/37)} {printf("%12s\t", d/37)}{printf("%12s\n", e/37) }' final.dat > freq.dat	
