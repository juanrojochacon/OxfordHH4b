#!/bin/sh

# takes n arguments, each a directory, the histogram dat files from which are to be passed to plotHist.
# one argument, the directory "analysis/sample" from which you want to plot the hists

if [ -z "$OUTPUTDIRLOC" ]; then
    echo "OUTPUTDIRLOC not set, defaulting to directory plots"
    export OUTPUTDIRLOC="plots"
fi

if [ ! -d $OUTPUTDIRLOC/ ]
then
mkdir $OUTPUTDIRLOC
fi

#export FILELIST=`for arg in $@; do echo $arg"/$dat";done`
#echo "Full list is " $FILELIST

for dat in `ls $1 | grep histo | grep dat | grep -v ptHptH | grep -v mHmH`
do {
export FILELIST=`for arg in $@; do echo $arg"/$dat";done`
#echo "Full list is " $FILELIST
python FIKRI_plotHist.py $FILELIST
mv histo.pdf $OUTPUTDIRLOC/`echo $dat | sed 's/.dat/.pdf/'`
#echo $dat
} done


# for dat in `ls $1 | grep histo | grep dat | grep ptHptH`
# do {
# for arg in $@
# do {
# #echo $arg"/$dat";done
# python plot2DHist.py $arg"/"$dat
# mv histo2D.pdf $OUTPUTDIRLOC/`echo $arg$dat | sed 's/\/h/_h/g' | sed 's/.dat/.pdf/'`
# } done
# #echo $dat
# } done



# for dat in `ls $1 | grep histo | grep dat | grep mHmH`
# do {
# for arg in $@
# do {
# #echo $arg"/$dat";done
# python plot2DHist.py $arg"/"$dat
# mv histo2D.pdf $OUTPUTDIRLOC/`echo $arg$dat | sed 's/\/h/_h/g' | sed 's/.dat/.pdf/'`
# } done
# #echo $dat
# } done

