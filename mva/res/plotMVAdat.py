#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
from math import sqrt
import sys, os, math
import numpy

# HL-LHC luminosity
hl_lhc_lumi=3000

# Setup ROC plot
roc, rocax = plt.subplots()
rocax.plot([0,1],[1,0], color='grey')

rocax.set_ylabel("Background rejection")
rocax.set_xlabel("Signal efficiency")

# Gridlines
rocax.xaxis.grid(True)
rocax.yaxis.grid(True)

# Setup s/b plot
sb, sbax = plt.subplots()
sbax.set_ylabel("S/B")
sbax.set_xlabel("NN Discriminant")

# Setup s/sqrt(b) plot
ssb, ssbax = plt.subplots()
ssbax.set_ylabel("S/sqrt(B)")
ssbax.set_xlabel("NN Discriminant")

# Setup N_evt plots
nev, nevax = plt.subplots()
nevax2 = nevax.twinx()
nevax.set_ylabel("Number of signal events")
nevax2.set_ylabel("Number of background events")
nevax.set_xlabel("NN Discriminant")

# Gridlines
sbax.xaxis.grid(True)
sbax.yaxis.grid(True)
ssbax.xaxis.grid(True)
ssbax.yaxis.grid(True)
nevax.xaxis.grid(True)
nevax.yaxis.grid(True)

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

	events = []	

	for line in infile:
		# Full event for s/b plot purposes
		events.append(line.split()[1:4])

		if line.split()[1] == '0':
			bkgprob.append( float(line.split()[3]) )
		else:
			sigprob.append( float(line.split()[3]) )

	bins = numpy.linspace(0, 1, 20)
	fig, ax = plt.subplots()

	ax.hist(bkgprob, bins, color='b', alpha=0.5, normed=True, label = "Background events")
	ax.hist(sigprob, bins, color='r', alpha=0.5, normed=True, label = "Signal events")
	ax.set_ylim([0,8])

	ax.set_xlabel("Neural network response")

	# Legend
	legend = ax.legend(fontsize=10, loc='best')
	legend.get_frame().set_alpha(0.7)

	numpoints = str( len(bkgprob) + len(sigprob) ) + " events: " + str(len(sigprob)) + " signal, " + str(len(bkgprob)) + " background."
	fig.text(0.13,0.92,numpoints, fontsize=12)
	fig.text(0.13,0.96,"MVA: "+ basename, fontsize=12)

	figname = basename + "_hist.pdf"
	fig.savefig(figname)

	#### ROC Curve and S/B plot

	thresholds = numpy.linspace(0, 1, 50)
	falsepos = []
	truepos = []

	soverb = []
	soversb = []

	nbkg = []
	nsig = []

	for th in thresholds:
		print "Calculating threshold: ", th
		
		# S/B signal and background weights
		signalwgt = 0
		backgdwgt = 0

		# ROC true and false positives
		tp = 0
		fp = 0

		for evt in events:
			signal = bool(int(evt[0]))
			weight = float(evt[1])
			discriminant = float(evt[2])

			if discriminant > th:
				if signal:
					signalwgt = signalwgt + weight # signal weight
					tp = tp+1 #ROC true positive
				else:
					backgdwgt = backgdwgt + weight # background weight
					fp = fp+1 # ROC false positive

		falsepos.append(1- fp/float(len(bkgprob)))
		truepos.append(tp/float(len(sigprob)))

		nsig.append(hl_lhc_lumi*signalwgt)
		nbkg.append(hl_lhc_lumi*backgdwgt)

		if backgdwgt == 0:
			soverb.append(0)
			soversb.append(0)
		else:
			soverb.append(signalwgt/backgdwgt)
			soversb.append((hl_lhc_lumi*signalwgt)/sqrt(hl_lhc_lumi*backgdwgt))

	rocax.plot(truepos, falsepos, label = basename)
	sbax.plot(thresholds, soverb, label = basename)
	ssbax.plot(thresholds, soversb, label = basename)

	nevax.plot(thresholds, nsig, color='r', label=basename+"_signal")
	nevax2.plot(thresholds, nbkg, color='b', label=basename+"_background")

# Legend
rlegend = rocax.legend(loc='best')
rlegend.get_frame().set_alpha(0.8)

# Legend
slegend = sbax.legend(loc='best')
slegend.get_frame().set_alpha(0.8)

# Legend
sslegend = ssbax.legend(loc='best')
sslegend.get_frame().set_alpha(0.8)

#NeV colors
for tl in nevax.get_yticklabels():
    tl.set_color('r')
for tl in nevax2.get_yticklabels():
    tl.set_color('b')

x1,x2,y1,y2 = ssbax.axis()
ssbax.axis((x1,x2,0,10))

roc.savefig('roc.pdf')
sb.savefig('sb.pdf')
ssb.savefig('ssb.pdf')
nev.savefig('nev.pdf')






