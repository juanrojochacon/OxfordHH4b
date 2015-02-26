#!/bin/sh

# takes n arguments, each a directory, the histogram dat files from which are to be passed to plotHist.

export OUTPUTDIRLOC="plots"

if [ ! -d $OUTPUTDIRLOC/ ]
then
mkdir $OUTPUTDIRLOC
fi

#export FILELIST=`for arg in $@; do echo $arg"/$dat";done`
#echo "Full list is " $FILELIST

for dat in `ls $1 | grep histo | head -2`
do {
export FILELIST=`for arg in $@; do echo $arg"/$dat";done`
echo "Full list is " $FILELIST
python plotHist.py $FILELIST
mv histo.pdf $OUTPUTDIRLOC/`echo $dat | sed 's/.dat/.pdf/'`
#echo $dat
} done

