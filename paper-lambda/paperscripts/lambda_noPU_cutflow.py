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

######################### Read in cutflows #########################

cutflows = {}
for lvalue in lvalues:
	folder="diHiggs_LAM"+str(lvalue)+"/"
	for regime in regimes:
		filename = root + folder + "histo_CF_"+regime+".yoda" 

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
			cutflows[folder+regime] = cutflow

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
			folder="diHiggs_LAM"+str(lvalue)+"/"
			xsec.append(cutflows[folder+regime][cut])
		ax.plot(lvalues, xsec, color=colours[icol])
		icol = icol + 1
	fig.savefig(regime+'.pdf')
	print(regime+'.pdf exported')

for regime in regimes:
	table = []
	for lvalue in lvalues:
		folder="diHiggs_LAM"+str(lvalue)+"/"
		table.append(cutflows[folder+regime])
	tab_trans = [list(i) for i in zip(*table)]
	tab = tabulate(tab_trans, lvalues,numalign="center", floatfmt=".2E")
	print("\n"+regime + " cutflow:")
	print(tab)


