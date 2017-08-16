#!/bin/zsh

mkdir -p ntuples
for file in  `find . -wholename "./SHERPA_QCD*/fullNTuple-*-tag.dat"`; do 
	ar=(${(s:/:)file})
	echo "Merging $ar..."
	cp $file ntuples/${ar[2]#SHERPA_QCD}-${ar[3]}
done
echo "Creating tarball"
tar acf ntuples.tar.xz ntuples

