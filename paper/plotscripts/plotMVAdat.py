#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
from math import sqrt
import sys, os, math
import numpy

print "Warning: This script might take a minute or two depending on the density of thresholds"

################################ Settings ###################################

# Source files
datafiles=[ "~/Dropbox/HH4bMC/mva/nn_13X5X3X1_500000-Gen_resNTuple.dat",
			"~/Dropbox/HH4bMC/mva/nn_13X5X3X1_500000-Gen_intNTuple.dat",
			"~/Dropbox/HH4bMC/mva/nn_17X5X3X1_500000-Gen_bstNTuple.dat"]

datanames=[ "Resolved",
			"Intermediate",
			"Boosted"]


# HL-LHC luminosity
hl_lhc_lumi=3000

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']

# Histogram names
Histout = "disc"  # Histograms for discriminant
ROCout = "roc"	  # ROC curve
SBout = "sb"	  # S/B curve
SSBout = "ssb"	  # S/Sqrt(B) curve
NeVout = "nev"	  # Number of events (twin axes)
NeV2out = "nev2"  # Number of events (single axis)

#############################################################################

# Plot the individual discriminant histograms
def plotDiscriminantHisto(name, signal, background):
	bins = numpy.linspace(0, 1, 20) # Binning density
	fig, ax = plt.subplots()

	ax.hist(signal, bins, color='r', alpha=0.5, normed=True, label = "Signal events")
	ax.hist(background, bins, color='b', alpha=0.5, normed=True, label = "Background events")

	ax.set_ylim([0,8])
	ax.set_xlabel("Neural network response")

	# Legend
	legend = ax.legend(fontsize=10, loc='best')
	legend.get_frame().set_alpha(0.7)

	numpoints = str( len(background) + len(signal) ) + " events: " + str(len(signal)) + " signal, " + str(len(background)) + " background."
	fig.text(0.13,0.92,numpoints, fontsize=12)
	fig.text(0.13,0.96,"MVA: "+ name, fontsize=12)

	figname = name + "_"+Histout+".pdf"
	fig.savefig(figname)

############################ Plot setup ########################################

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

nev3, nevax3 = plt.subplots()
nevax3.set_ylabel("Number of events")
nevax3.set_xlabel("NN Discriminant")

# Gridlines
sbax.xaxis.grid(True)
sbax.yaxis.grid(True)
ssbax.xaxis.grid(True)
ssbax.yaxis.grid(True)
nevax.xaxis.grid(True)
nevax.yaxis.grid(True)
nevax3.yaxis.grid(True)
nevax3.yaxis.grid(True)

nevax3.set_yscale('log')

######################## Reading data ##############################

for idat in xrange(0,len(datafiles)):
	infilenm = os.path.expanduser(datafiles[idat])
	basename = datanames[idat]

	# Verify paths
	if os.path.exists(infilenm) == False:
	  	print "Error: source file" + infilenm + " not found!"
	  	sys.exit()

	infile = open(infilenm, 'rb')
	print "Processing " + infilenm + " ..."


	# Setup arrays
	bkgprob = [] # Background point discriminant
	sigprob = [] # Signal point discriminant
	events = []	 # Full event info

	# Read data
	for line in infile:
		events.append(line.split()[1:4])
		if line.split()[1] == '0':
			bkgprob.append( float(line.split()[3]) )
		else:
			sigprob.append( float(line.split()[3]) )

	#### ROC Curve and S/B plot
	thresholds = numpy.linspace(0, 1, 5)
	falsepos = []
	truepos = []

	soverb = []
	soversb = []

	nbkg = []
	nsig = []

	for th in thresholds:
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

	# Plot discriminant histogram
	plotDiscriminantHisto(basename, sigprob, bkgprob)

	# Plot ROC curve, s/b, s/sqrt(b)
	rocax.plot(truepos, falsepos, color=colours[idat], label = basename)
	sbax.plot(thresholds, soverb, color=colours[idat], label = basename)
	ssbax.plot(thresholds, soversb, color=colours[idat], label = basename)

	# Plot number of events per discriminant cut
	nevax.plot(thresholds, nsig, color=colours[idat], linestyle='--', label=basename)
	nevax2.plot(thresholds, nbkg, color=colours[idat], label=basename)

	nevax3.plot(thresholds, nsig, color=colours[idat], linestyle='--')
	nevax3.plot(thresholds, nbkg, color=colours[idat], label=basename)

################################### Finish up ########################################

# Legends
rlegend = rocax.legend(loc='best')
rlegend.get_frame().set_alpha(0.8)
slegend = sbax.legend(loc='best')
slegend.get_frame().set_alpha(0.8)
sslegend = ssbax.legend(loc='best')
sslegend.get_frame().set_alpha(0.8)
nevlegend = nevax2.legend(loc='best')
nevlegend.get_frame().set_alpha(0.8)
nev3legend = nevax3.legend(loc='best')
nev3legend.get_frame().set_alpha(0.8)

x1,x2,y1,y2 = ssbax.axis()
ssbax.axis((x1,x2,0,4))

roc.savefig(ROCout+'.pdf')
sb.savefig(SBout+'.pdf')
ssb.savefig(SSBout+'.pdf')
nev.savefig(NeVout+'.pdf')
nev3.savefig(NeV2out+'.pdf')






