#! /usr/local/Cellar/python/2.7.10/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
from math import sqrt
import sys, os, math
import numpy

print "Warning: This script might take a minute or two depending on the density of thresholds"

################################ Settings ###################################

# Source files
datafiles=[ "../plotdata/results_SK_PU40/MVA/nn_13X5X3X1_50000-Gen_resNTuple.dat",
            "../plotdata/results_SK_PU40/MVA/nn_17X5X3X1_50000-Gen_intNTuple.dat",
            "../plotdata/results_SK_PU40/MVA/nn_21X5X3X1_50000-Gen_bstNTuple.dat"]

datanames=[ "Resolved","Intermediate","Boosted"]
titlenames=[ r"Resolved Category,  $ n_{PU} =40$+SK",
             r"Intermediate Category,  $ n_{PU} =40$+SK",\
             r"Boosted Category,  $ n_{PU}=40$+SK"]

datanamessignal=[ "Signal Res","Signal Int","Signal Boost"]
datanamesback=[ "Background Res","Background Int","Background Boost"]


# HL-LHC luminosity
hl_lhc_lumi=3000

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']

linestyles = ['dashed','dotted','solid']

# Histogram names
Histout = "disc_SKPU40"  # Histograms for discriminant
ROCout = "roc_SKPU40"	  # ROC curve
SBout = "sb_SKPU40"	  # S/B curve
SSBout = "ssb_SKPU40"	  # S/Sqrt(B) curve
NeVout = "nev_SKPU40"	  # Number of events (twin axes)
NeV2out = "nev2_SKPU40"  # Number of events (single axis)

# Rescaling factor for signal
resc = 0.41

 
#############################################################################

# Plot the individual discriminant histograms
def plotDiscriminantHisto(name, signal, background,title):
        plt.rcParams.update({'font.size': 16})
	bins = numpy.linspace(0, 1, 20) # Binning density
	fig, ax = plt.subplots()

	ax.hist(signal, bins, color='r', alpha=0.8, normed=True, label = "Signal",linewidth=2)
	ax.hist(background, bins, color='b',fill=False, alpha=0.8, edgecolor='b',normed=True, label = "Background",hatch="//",linewidth=2)

	ax.set_ylim([0,10])
	ax.set_xlabel("ANN Output")
        ax.set_ylabel("a. u.")

	# Legend
	legend = ax.legend(fontsize=18, loc='best')
	legend.get_frame().set_alpha(0.8)

	numpoints = str( len(background) + len(signal) ) + " events: " + str(len(signal)) + " signal, " + str(len(background)) + " background."
	fig.text(0.23,0.93,title, fontsize=19)

	figname = name + "_"+Histout+".pdf"
	fig.savefig(figname)

############################ Plot setup ########################################

# Setup ROC plot
roc, rocax = plt.subplots()
rocax.plot([0,1],[1,0], color='grey')

rocax.set_ylabel("Background rejection rate",fontsize=17)
rocax.set_xlabel("Signal efficiency",fontsize=17)

# Gridlines

# Setup s/b plot
sb, sbax = plt.subplots()
sbax.set_ylabel("$S/B$",fontsize=19)
sbax.set_xlabel("ANN output cut")
sbax.set_ylim([0.0001,1])
sbax.set_xlim([0,0.9])

# Setup s/sqrt(b) plot
ssb, ssbax = plt.subplots()
ssbax.set_ylabel("$S/\sqrt{B}$",fontsize=19)
ssbax.set_xlabel("ANN output cut")
ssbax.set_xlim([0,0.9])

# Setup N_evt plots
nev, nevax = plt.subplots()
nevax2 = nevax.twinx()
nevax.set_ylabel("Number of signal events")
nevax2.set_ylabel("Number of background events")
nevax.set_xlabel("ANN output cut")

nev3, nevax3 = plt.subplots()
nevax3.set_ylabel("$N_{ev}$ at HL-LHC after MVA cut")
nevax3.set_xlabel("ANN output cut")
nevax3.set_xlim([0,0.9])

sbax.set_yscale('log')
nevax3.set_yscale('log')
nevax3.set_ylim([1,1e9])

######################## Reading data ##############################

for idat in xrange(0,len(datafiles)):
	infilenm = os.path.expanduser(datafiles[idat])
	basename = datanames[idat]
        titlename = titlenames[idat]
        basenamesignal = datanamessignal[idat]
        basenameback = datanamesback[idat]

	# Verify paths
	if os.path.exists(infilenm) == False:
	  	print "Error: source file" + infilenm + " not found!"
	  	sys.exit()

	infile = open(infilenm, 'rb')
	print "Processing " + infilenm + " ..."


	# Setup arrays
	bkgprob = [] # Background point discriminant
	sigprob = [] # Signal point discriminant
	total_bkgweight = 0 # Total background weight
	total_sigweight = 0 # Total signal weight
	events = []	 # Full event info

	# Read data
	for line in infile:
		events.append(line.split()[1:4])
		if line.split()[1] == '0':
			bkgprob.append( float(line.split()[3]) ) # Discriminant
			total_bkgweight = total_bkgweight + float(line.split()[2]) # Weight
		else:
			sigprob.append( float(line.split()[3]) ) # Discriminant
			total_sigweight = total_sigweight + float(line.split()[2]) # Weight

        #####################################################                
        # Apply overall rescaling of the signal weight
        total_sigweight = resc * total_sigweight
        #####################################################                

	#### ROC Curve and S/B plot
	thresholds = numpy.linspace(0, 1, 40)
	falsepos = []
	truepos = []

	soverb = []
	soversb = []

	nbkg = []
	nsig = []

        print "th,   s_over_b,      s_over_sb,    Nev_sig,    Nev_back"
	for th in thresholds:
		# S/B signal and background weights
		signalwgt = 0
		backgdwgt = 0

		for evt in events:
			signal = bool(int(evt[0]))
			weight = float(evt[1])
			discriminant = float(evt[2])

			if discriminant > th:
				if signal:
					signalwgt = signalwgt + weight # signal weight
				else:
					backgdwgt = backgdwgt + weight # background weight

                                     
                #####################################################                
                # Apply overall rescaling of the signal weight                        
                signalwgt = resc * signalwgt
                #####################################################                
		falsepos.append(1 - backgdwgt/total_bkgweight)
		truepos.append(signalwgt/total_sigweight)

		nsig.append(hl_lhc_lumi*signalwgt)
		nbkg.append(hl_lhc_lumi*backgdwgt)

		if backgdwgt == 0:
			soverb.append(0)
			soversb.append(0)
		else:
			soverb.append(signalwgt/backgdwgt)
			soversb.append((hl_lhc_lumi*signalwgt)/sqrt(hl_lhc_lumi*backgdwgt))
                        s_over_b = signalwgt/backgdwgt
                        s_over_sb = (hl_lhc_lumi*signalwgt)/sqrt(hl_lhc_lumi*backgdwgt)
                        print th, s_over_b," ", s_over_sb, " ",hl_lhc_lumi*signalwgt, " ",hl_lhc_lumi*backgdwgt

	# Plot discriminant histogram
	plotDiscriminantHisto(basename, sigprob, bkgprob,titlename)

	# Plot ROC curve, s/b, s/sqrt(b)
	rocax.plot(truepos, falsepos, color=colours[idat],linestyle=linestyles[idat], label = basename,linewidth=2.4)
	sbax.plot(thresholds, soverb, color=colours[idat], linestyle=linestyles[idat], label = basename,linewidth=2.4)
	ssbax.plot(thresholds, soversb, color=colours[idat],linestyle=linestyles[idat], label = basename,linewidth=2.4)

	# Plot number of events per discriminant cut
	nevax.plot(thresholds, nsig, color=colours[idat], linestyle='--', label=basename,linewidth=2.4)
	nevax2.plot(thresholds, nbkg, color=colours[idat], label=basename,linewidth=2.4)

        nevax3.plot(thresholds, nsig, color=colours[idat], label=basenamesignal, linestyle='--',linewidth=2.4)
	nevax3.plot(thresholds, nbkg, color=colours[idat], label=basenameback,linewidth=2.4)
        

################################### Finish up ########################################

# Legends
plt.rcParams.update({'font.size': 19})
rlegend = rocax.legend(loc='best')
rlegend.get_frame().set_alpha(0.8)
slegend = sbax.legend(loc='best')
slegend.get_frame().set_alpha(0.8)
sslegend = ssbax.legend(loc='best')
sslegend.get_frame().set_alpha(0.8)
nevlegend = nevax2.legend(loc='best')
nevlegend.get_frame().set_alpha(0.8)
nev3legend = nevax3.legend(loc='best',fontsize=13)
nev3legend.get_frame().set_alpha(0.8)

x1,x2,y1,y2 = ssbax.axis()
ssbax.axis((x1,x2,0,4.2))

roc.text(0.29,0.94,r"HL-LHC, PU40+SK+Trim", fontsize=20)
roc.savefig(ROCout+'.pdf')

sb.text(0.29,0.94,r"HL-LHC, PU40+SK+Trim", fontsize=20)
sb.savefig(SBout+'.pdf')

ssb.text(0.29,0.94,r"HL-LHC, PU40+SK+Trim", fontsize=20)
ssb.savefig(SSBout+'.pdf')

nev3.text(0.29,0.94,r"HL-LHC, PU40+SK+Trim", fontsize=20)
nev3.savefig(NeV2out+'.pdf')






