#! /usr/bin/python
import sys, os, math

desc = []
sumvals = []
for idat in xrange(1,len(sys.argv)):

	infilenm = sys.argv[idat]
	basename = os.path.splitext(infilenm)[0]

	# Verify paths
	if os.path.exists(infilenm) == False:
	  	print "Error: source file" + infilenm + " not found!"
	  	sys.exit()
	infile = open(infilenm, 'rb')

	icounter = 0
	for line in infile:
		meta = " ".join(line.split()[0:3])
		if len(desc) <= icounter:
			desc.append(meta)
			sumvals.append(float(line.split()[3]))
		else:
			if desc[icounter] != meta:
				print("error " + desc[icounter] +" "+ meta)
			else:
				sumvals[icounter] = sumvals[icounter] + float(line.split()[3])
		icounter = icounter + 1


for iline in xrange(0,len(desc)):
	print(desc[iline] + " "+str(sumvals[iline]/(len(sys.argv)-1)))
