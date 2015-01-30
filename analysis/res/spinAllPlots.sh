#!/bin/sh

# takes one argument, the directory "analysis/sample" from which you want to plot the hists

export OUTPUTDIRLOC="plots"

if [ ! -d $OUTPUTDIRLOC/ ]
then
mkdir $OUTPUTDIRLOC
fi

for dat in `ls $1 | grep histo`
do {
python plotHist.py $1/$dat
mv histo.pdf $OUTPUTDIRLOC/`echo $dat | sed 's/.dat/.pdf/'`
#echo $dat
} done

