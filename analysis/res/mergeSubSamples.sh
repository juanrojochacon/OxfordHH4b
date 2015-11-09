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
    	COUNT=0
        for f in $d*NTuple.*.dat; do
			filename=$(basename $f)
			basename=${filename:0:9}
        	echo "NTuple Merging " $f " > " $basename.dat

			if [[ $COUNT == 0 ]]; then
				cat $f > $d$basename.dat
			else
				cat $f | awk ' NR>1 {print;}' >> $d$basename.dat
			fi

        	COUNT=$((COUNT+1))
        done
        
        # Cleanup
        rm $d*NTuple.*.dat

		# Now merge histograms into total background 
	    if [[ $d != *diHiggs* ]]; then
			for f in $d*.yoda; do
				if [[ $SCOUNT == 1 ]]; then
					echo "Total Background Generation " $f
					cp $f $1background/$( basename $f)
				else
					echo "Total Background Merging " $f
					yodamerge -o $1background/$( basename $f) $f $1background/$( basename $f)
				fi
			done

			SCOUNT=$(( $SCOUNT+1 ))
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

# Now merge background and signal nTuples into totals
# Start by copying diHiggs nTuples into total
for f in $1diHiggs/*NTuple.dat
do
	cp $f $1
done

# Now merge in all background samples
for d in $1/*/
do
  	if [[ $d != *diHiggs* ]]; then
                if [[ $d != *background* ]]; then
                        for f in $d*NTuple.dat; do
                        	echo "Total Ntuple merging " $f
							filename=$(basename $f)
                        	cat $f | awk ' NR>1 {print;}' >> $1$filename 
                        done
                fi
        fi
done


