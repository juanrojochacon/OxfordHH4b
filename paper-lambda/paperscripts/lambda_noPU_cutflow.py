#!/usr/bin/python

import math
import yoda
from tabulate import tabulate
from matplotlib import pyplot as plt

##################################################################

root = "../plotdata/noPU_results/"

lvalues = [ -7, -5, -3, -1, 0, 1, 2, 3, 5, 7, 9 ]
regimes = ['res', 'inter','boost']

Cnum = ['C0', 'C1a', 'C1b', 'C1c', 'C1d', 'C1e', 'C2' ]

lumi = 3000.0
sysErr = 0.1 # 10%

def fname(lval):
	return "diHiggs_LAM"+str(lval)+"/"

######################### Read in cutflows #########################

cutflows = {}
for lvalue in lvalues:
	for regime in regimes:
		filename = root + fname(lvalue) + "histo_CF_"+regime+".yoda" 

		## Load and sort data objects
		aos = yoda.read(filename)
		for aopath, ao in aos.iteritems():
			if type(ao) != yoda.Histo1D:
				print "Error: AnalysisObject is not a Histo1D"
				quit();

			# Add cutflow
			cutflow = []
			for bin in ao.bins:
				cutflow.append(bin.sumW)
			cutflows[fname(lvalue)+regime] = cutflow

######################### Make plots #########################

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
for regime in regimes:
	fig, ax = plt.subplots()
	ax.set_ylabel("HH xSec")
	ax.set_xlabel("lambda")
	ax.set_yscale('log')
	icol=0
	#axes.set_xlim([xmin,xmax])
	fig.suptitle(regime)
	ax.set_ylim([1E-2,1E3])
	for cut in range(0,len(Cnum)):
		xsec = []
		for lvalue in lvalues:
			xsec.append(cutflows[fname(lvalue)+regime][cut])
		ax.plot(lvalues, xsec, color=colours[icol])
		icol = icol + 1
	fig.savefig(regime+'.pdf')
	print(regime+'.pdf exported')

chi2tab_all = [] # All chi2 tables
for regime in regimes:
	table = []
	for lvalue in lvalues:
		table.append(cutflows[fname(lvalue)+regime])

	# Table for printing
	tab_trans = [list(i) for i in zip(*table)]
	tab = tabulate(tab_trans, lvalues,numalign="center", floatfmt=".2E")
	print("\n"+regime + " cutflow:")
	print(tab)

	# Now chi2
	chi2table = []
	for lvalue in lvalues:
		chi2vals = []
		for cut in range(0,len(Cnum)):
			SigLam = cutflows[fname(lvalue)+regime][cut] # Cross-section at lambda
			SigSM = cutflows[fname(1)+regime][cut]		 # Cross-section at lambda=1

			if SigLam != 0: # Intermediate has a unfilled cut
				err1 = math.sqrt(SigSM/lumi)
				err2 = math.sqrt(SigLam/lumi)
				chi2 = pow(SigSM-SigLam,2)/(pow(err1,2) +pow(err2,2)+ pow(sysErr*SigSM,2))
				chi2vals.append(chi2)
			else: 
				chi2vals.append(-1)
		chi2table.append(chi2vals)

	chi2tab_trans = [list(i) for i in zip(*chi2table)]
	chi2tab_all.append(chi2tab_trans[len(Cnum)-1]) # Append final cutflow 

# Prepare chi2 plots
fig, ax = plt.subplots()
ax.set_ylabel("$\chi^2$")
ax.set_xlabel("$\lambda$")
ax.set_ylim([0,10])
fig.suptitle("$\chi^2$ profile for all topologies L="+str(lumi)+"fb$^{-1}$")

# Print out final values
ireg=0
labelname=["Resolved", "Intermediate", "Boosted"]
for chi2plot in chi2tab_all:
	ax.plot(lvalues, chi2plot, label=labelname[ireg])
	ireg=ireg+1

# Now add the legend with some customizations.
legend = ax.legend(loc='best', shadow=True)
fig.savefig('chi2.pdf')
