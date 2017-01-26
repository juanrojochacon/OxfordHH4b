#!/bin/bash
# This script merges all the samples in a specified folder
# into combined background samples and NTuples.
rm -rf $1/background
mkdir $1/background

for d in $1*/
do
	if [[ $d != *background* ]]; then
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

# Convert yoda files in background to flat
for f in $1background/*.yoda
do
	filename=$( basename $f )
	filename="${filename%.*}".dat
	yoda2flat $f $1background/$filename
done
