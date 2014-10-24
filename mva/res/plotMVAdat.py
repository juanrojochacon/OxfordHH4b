#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy


# Setup ROC plot
roc, rocax = plt.subplots()
rocax.plot([0,1],[0,1], color='grey')

rocax.set_xlabel("False positive rate")
rocax.set_ylabel("True positive rate")

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

	for line in infile:
		if line.split()[1] == '0':
			bkgprob.append( float(line.split()[2]) )
		else:
			sigprob.append( float(line.split()[2]) )

	bins = numpy.linspace(0, 1, 20)
	fig, ax = plt.subplots()

	ax.hist(bkgprob, bins, color='b', alpha=0.5, normed=True, label = "Background events")
	ax.hist(sigprob, bins, color='r', alpha=0.5, normed=True, label = "Signal events")
	ax.set_ylim([0,8])


	legend = ax.legend(loc='best')
	legend.get_frame().set_alpha(0.8)

	ax.set_xlabel("Neural network response")

	numpoints = str( len(bkgprob) + len(sigprob) ) + " events: " + str(len(sigprob)) + " signal, " + str(len(bkgprob)) + " background."
	fig.text(0.13,0.92,numpoints, fontsize=12)
	fig.text(0.13,0.96,"MVA: "+ basename, fontsize=12)

	figname = basename + "_hist.pdf"
	fig.savefig(figname)

	#### ROC Curve

	thresholds = numpy.linspace(0, 1, 200)
	falsepos = []
	truepos = []

	for th in thresholds:
		fp = 0
		for bkg in bkgprob:
			if bkg > th:
				fp = fp+1
		falsepos.append(fp/float(len(bkgprob)))

		tp = 0
		for sig in sigprob:
			if sig > th:
				tp = tp+1
		truepos.append(tp/float(len(sigprob)))


	rocax.plot(falsepos,truepos, label = basename)

# Gridlines
rocax.xaxis.grid(True)
rocax.yaxis.grid(True)

# Legend
legend = rocax.legend(loc='best')
legend.get_frame().set_alpha(0.8)
roc.savefig('roc.pdf')




