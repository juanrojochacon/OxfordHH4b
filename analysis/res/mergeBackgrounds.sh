#!/bin/bash

rm -rf $1/background
mkdir $1/background

SCOUNT=1

for d in $1/*/
do
	if [[ $d != *diHiggs* ]]; then
		if [[ $d != *background* ]]; then
			for f in $d*.yoda; do
				if [[ $SCOUNT == 1 ]]; then
					echo "Copying " $f
					cp $f $1/background/$( basename $f)
				else
					echo "Merging " $f
					yodamerge -o $1/background/$( basename $f) $f $1/background/$( basename $f)
				fi
			done

			SCOUNT=$(( $SCOUNT+1 ))
		fi
	fi
done
