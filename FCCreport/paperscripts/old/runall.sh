#!/bin/bash

for f in *.py; do 
	echo "Processing $f.."; 
	python $f
done

