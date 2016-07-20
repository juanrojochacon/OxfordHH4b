#!/bin/bash
cd ../../mva/
make clean && make
cd ../paper-lambda/plotdata
cp ../../mva/apply_mva ./

SOURCE='4b'
ANNCUT=$1

INPUTDIR="./baseline_noPU_sideband_14/"
TARGETDIR="./"$SOURCE"_res/"
rm -rf $TARGETDIR
mkdir $TARGETDIR

 LVALS=(-7 -5 -3 -1 0 0.2 0.4 0.6 0.8 1 1.2 1.4 1.8 2 3 5 7 9)
 NTUPLES=( 'bst' 'int' 'res' )
 NTUPLES_LONG=( 'boost' 'inter' 'res' )
 for LAMBDA in "${LVALS[@]}"; do
 	LFILE=${LAMBDA/./_}
 	DIR=$INPUTDIR"/diHiggs_LAM$LFILE/"
 	for i in $(seq 0 2); do
 		NTUP=${NTUPLES[$i]}
 		NTUPL=${NTUPLES_LONG[$i]}
 		# Setup cutflow
		echo $TARGETDIR$SOURCE"_"$LFILE"_"$NTUP".out"
 		sed -n '6,8p' $INPUTDIR"diHiggs_LAM"$LFILE"/histo_CF_"$NTUPL".dat" > $TARGETDIR$SOURCE"_"$LFILE"_"$NTUP".out"
 		echo $DIR$NTUP"NTuple.dat" $INPUTDIR$NTUP".par"
 		./apply_mva $DIR$NTUP"NTuple.dat" $INPUTDIR$NTUP".par" $ANNCUT 
 		cat $TARGETDIR$SOURCE"_"$LFILE"_"$NTUP".out" ./results.dat > $TARGETDIR$SOURCE"_"$LFILE"_"$NTUP".dat"
 		rm results.dat
 	done
 done

rm $TARGETDIR*.out
rm apply_mva
