#! /bin/bash
for f in *.dat
do
	./plotHist.py ${f}
 #yoda2flat ${f}
 #yodaplot ${f%%.*}".dat"
done


