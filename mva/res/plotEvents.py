#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
from math import sqrt
import sys, os, math
import numpy

# HL-LHC luminosity
hl_lhc_lumi=3000
bins = numpy.linspace(0, 1, 20)
for idat in xrange(1,len(sys.argv)):

	infilenm = sys.argv[idat]
	basename = os.path.splitext(infilenm)[0]

	# Verify paths
	if os.path.exists(infilenm) == False:
	  	print "Error: source file" + infilenm + " not found!"
	  	sys.exit()

	infile = open(infilenm, 'rb')
	print "Processing " + infilenm + " ..."

	bkgprob = []
	sigprob = []
	bkgwgts = []
	sigwgts = []

	for line in infile:
		if line.split()[1] == '0':
			bkgprob.append( float(line.split()[3]) )
			bkgwgts.append( hl_lhc_lumi*float(line.split()[2]) )
		else:
			sigprob.append( float(line.split()[3]) )
			sigwgts.append( hl_lhc_lumi*float(line.split()[2]) )

	fig, ax = plt.subplots()

	totprob = bkgprob + sigprob
	totwgts = bkgwgts + sigwgts

	ax.hist(totprob, bins, weights = totwgts, color='r', alpha=1, label = "Signal events")
	ax.hist(bkgprob, bins, weights = bkgwgts, color='b', alpha=1, label = "Background events")
	#ax.hist(sigprob, bins, weights = sigwgts, color='r', alpha=1, label = "Signal events")
	ax.set_ylim([1E3,1E5])

	ax.set_xlabel("Neural network response")
	ax.set_yscale('symlog')

	# Legend
	legend = ax.legend(fontsize=10, loc='best')
	legend.get_frame().set_alpha(0.7)

	numpoints = str( len(bkgprob) + len(sigprob) ) + " events: " + str(len(sigprob)) + " signal, " + str(len(bkgprob)) + " background."
	fig.text(0.13,0.92,numpoints, fontsize=12)
	fig.text(0.13,0.96,"MVA: "+ basename, fontsize=12)

	figname = basename + "_hist.pdf"
	fig.savefig(figname)


