#! /usr/local/Cellar/python/2.7.10/bin/python
import math
import yoda
import os
import string

##################################################################

# input nTuple and NN discriminant

# Boosted category, PU80+SK+Trim
#infile = "~/Dropbox/HH4bMC/mva/results_SK_PU80_trim/bst_SKPU80_trim.dat"
# nnfile = "../plotdata/results_SK_PU80/MVA/nn_21X5X3X1_50000-Gen_SKPU80_bst.dat" 

# Resolved category, no PU
infile = "~/Dropbox/HH4bMC/mva/results_noPU/resNTuple.dat"
nnfile="../plotdata/results_noPU/MVA/nn_13X5X3X1_50000-Gen_noPU_res.dat"

infile = os.path.expanduser(infile)

# NN discriminant cut
ycut = 0.0
# ycut = 0.60 # optimal ANN cut in the resolved category wo PU

signal = 1 # switch for signal recognition (1 for signal, 0 for background)

######################### Read nTuple kinematic limits #########################

nTuple = open(infile, 'r')
kinvar = nTuple.readline().split()[4:] # List of kinematic variables
nkin = len(kinvar)

kinmin = [ float("inf") for i in kinvar ]
kinmax = [ 0 for i in kinvar ]

for line in nTuple:
	var = line.split()[3:]
	for i in xrange(0,nkin):
		kinmin[i] = min(kinmin[i], float(var[i]))
		kinmax[i] = max(kinmax[i], float(var[i]))

################################ Init Histograms ################################
histos = [ yoda.Histo1D(20, kinmin[i], kinmax[i], title=kinvar[i]) for i in xrange(0,nkin) ]
	
# Reopen nTuple
nTuple = open(infile, 'r')
nTuple.readline()

# Open NN discriminant file
NNout = open(nnfile, 'r')

for line in nTuple:
	var = line.split()[3:]
	disc = NNout.readline().split()

	assert line.split()[2] == disc[2] # weight check
	ynn =float(disc[3]) # NN discriminant

	if (ynn > ycut) and disc[1] == str(signal):
		for i in xrange(0,nkin):
			histos[i].fill(float(var[i]), float(disc[2]))

ext = '_y'+string.replace(str(ycut),'.','')+'_sig'+str(signal)+'.dat'
for i in xrange(0,nkin):
	yoda.core.writeFLAT(histos[i],'histo_'+kinvar[i]+ext)

