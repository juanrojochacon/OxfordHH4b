#!/bin/bash

rm -rf $1/background
mkdir $1/background

for d in $1/*/
do
	echo "Merging " $d
	for f in $d*.0.yoda; do
		FILENAME=$(basename $f)
		BASE="${FILENAME%.*}"
                BASE="${BASE%.*}"

                echo "Merging " $BASE

		yodamerge -o $d$BASE.dat $d$BASE.*.yoda
		rm $d$BASE.*.yoda
	done
done
