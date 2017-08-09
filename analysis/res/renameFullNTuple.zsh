#!/bin/zsh

mkdir -p ntuple
for file in  `find . -wholename "./SHERPA_QCD*/fullNTuple-*-tag.dat"`; do 
	ar=(${(s:/:)file})
	cp $file ntuple/${ar[2]#SHERPA_QCD}-${ar[3]}
done

