#!/bin/bash

for d in $1*/
do
	for f in $d*.2.yoda; do
		FILENAME=$(basename $f)
		BASE="${FILENAME%.*}"
        BASE="${BASE%.*}"

        echo "Merging " $d$BASE

		yodamerge -o $d$BASE.mrg $d$BASE.*.yoda
		rm $d$BASE.*.yoda
	done

	for f in $d*.mrg; do
		echo "moving "$f
                FILENAME=$(basename $f)
                BASE="${FILENAME%.*}"
		mv $f ./$d$BASE.yoda
		yoda2flat ./$d$BASE.yoda ./$d$BASE.dat
	done

	#Merge Ntuples
	if [[ $d != *background* ]]; then
    	COUNT=0
            for f in $d*NTuple.*.dat; do
            	echo "merging " $f 
				filename=$(basename $f)
				basename=${filename:0:9}

				if [[ $COUNT == 0 ]]; then
					cat $f > $d$basename.dat
				else
					cat $f | awk ' NR>1 {print;}' >> $d$basename.dat
				fi

            	COUNT=$((COUNT+1))
            done
    fi


done
