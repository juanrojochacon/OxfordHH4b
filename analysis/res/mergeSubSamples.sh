#!/bin/bash

rm -rf $1/background
mkdir $1/background

# First merge the various subsamples

SCOUNT=1
for d in $1*/
do
	if [[ $d != *background* ]]; then

		# Merge subsample yoda files
		for f in $d*.2.yoda; do
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

		# Now merge histograms into total background 
	    if [[ $d != *diHiggs* ]]; then
			for f in $d*.yoda; do
				TARGET=$1background/$( basename $f )
				if [ -f $TARGET ]; 
				then
					echo "Total Background Merging " $f
					yodamerge -o $TARGET $f $TARGET
				else
					echo "Total Background Generation " $f
					cp $f $TARGET
				fi
			done
		fi
    fi
done


# Convert yoda files in background to flat
for f in $1background/*.yoda
do
	filename=$( basename $f )
	filename="${filename%.*}".dat
	yoda2flat $f $1background/$filename
done


# Now merge totether nTuples
for d in $1/*/
do
    if [[ $d != *background* ]]; then
        for f in $d*NTuple.dat; do
			filename=$(basename $f)

        	if [ -f $1$filename ]; 
			then
            	echo "Total Ntuple merging " $f
            	cat $f | awk ' NR>1 {print;}' >> $1$filename 
            else
            	echo "Total NTuple generation " $f
            	cp $f $1$filename
        	fi
        done
    fi
done


