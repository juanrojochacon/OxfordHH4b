#!/usr/bin/env zsh
# This script merges subsample runs of all the samples
# in a specified directory

# Switch to using an alias instead of passing location of binary on the command line
alias yoda-merge=/home/stanislaus/yoda-merge
for d in $1*
do
    if [[ $d != *background* ]]; then
       #Merge Ntuples from Subsamples
        for f in $d/fullNTuple-<1-5>-tag.<0->.dat; do
            basename=${f:r:r}
            echo "NTuple Merging " $f " > " $basename.dat

            if [ -f $basename.dat ];
            then
		# remove first line (column headings)
                tail -n +2 $f >> $basename.dat
            else
                echo "Generating file"
                cp $f $basename.dat

            fi

        done

        # Cleanup
#        rm $d/fullNTuple-<1-5>-tag.*.dat
    fi
done
