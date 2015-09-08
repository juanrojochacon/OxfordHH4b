#Merge Ntuples
for f in $1/diHiggs/*NTuple.dat
do
	cp $f $1/
done

for d in $1/*/
do
  	if [[ $d != *diHiggs* ]]; then
                if [[ $d != *background* ]]; then
                        for f in $d*NTuple.dat; do
				filename=$(basename $f)
                        	cat $f | awk ' NR>1 {print;}' >> $1/$filename 
                        done
                fi
        fi
done

#cat $1/diHiggs/ntuple.res 






