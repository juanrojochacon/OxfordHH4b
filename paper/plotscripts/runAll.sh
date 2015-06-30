#!/bin/bash

for f in *.py; do 
 echo "Processing $f script.."; 
 python $f
done

mv ./*.pdf ../autoplots/
