#!/bin/bash
# This script merges subsample runs of all the samples
# in a specified directory

for d in $1*/
do
	if [[ $d != *background* ]]; then

		# Merge subsample yoda files
		for f in $d*.0.yoda; do
			FILENAME=$(basename $f)
			BASE="${FILENAME%.*}"
		        BASE="${BASE%.*}"

		        echo "Subsample Merging " $d$BASE
			yodamerge -o $d$BASE.mrg $d$BASE.*.yoda
			rm $d$BASE.*.yoda
		done

		# Sort out filenames for merged subsamples - process to .dat for plotHist
		for f in $d*.mrg; do
			echo "Merged Subsample Moving "$f
	        FILENAME=$(basename $f)
	        BASE="${FILENAME%.*}"
			mv $f ./$d$BASE.yoda
			yoda2flat ./$d$BASE.yoda ./$d$BASE.dat
		done

		#Merge Ntuples from Subsamples
        for f in $d*NTuple.*.dat; do
			filename=$(basename $f)
			basename=${filename:0:9}
        	echo "NTuple Merging " $f " > " $basename.dat

			if [ -f $d$basename.dat ];
			then
				cat $f | awk ' NR>1 {print;}' >> $d$basename.dat
			else
				echo "Generating file"
				cat $f > $d$basename.dat
			fi

        done
        
        # Cleanup
        rm $d*NTuple.*.dat
    fi
done
