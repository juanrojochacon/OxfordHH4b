#! /usr/local/Cellar/python/2.7.10/bin/python

import math
import yoda
from tabulate import tabulate

##################################################################


root = "../../paper/plotdata/results_noPU/"
labeltag="noPU"

#root = "../../paper/plotdata/results_SK_PU80/"
#labeltag="SKPU80"

folders = ['diHiggs/', 'background/', 'SHERPA_QCD4b/', 'SHERPA_QCD2b2j/', 'SHERPA_QCD4j/', 'SHERPA_QCDttbar/']
names = ['signal', 'total bkg', '4b', '2b2j', '4j', 'ttbar']

regimes = ['res', 'inter','boost']

Cnum = ['C0', 'C1a', 'C1b', 'C1c', 'C1d', 'C1e', 'C2' ]

lumi = 3000.0 # used for S/sqrt(B) only.
tablefmt = 'latex'

######################### Cutflow functions ########################

def Null(S,cutflow):
	return cutflow

def Nevt(S,cutflow):
	events = []
	for cut in cutflow:
		events.append(lumi*cut)
	return events

def SoverB(S,B):
	SoB = []
	for i in range(0,len(S)):
		if B[i] != 0:
			SoB.append(S[i]/B[i])
		else:
			SoB.append(0)
	return SoB

def SoverSqrtB(S,B):
	SosB = []
	for i in range(0,len(S)):
		if B[i] != 0:
			SosB.append(lumi*S[i]/math.sqrt(lumi*B[i]))
		else:
			SosB.append(0)
	return SosB

operations = {	"": Null,
				" Nevt": Nevt, 
				" S/B": SoverB, 
				" S/sqrt(B)": SoverSqrtB
			 }

######################### Read in cutflows #########################

cutflows = {}
for folder in folders:
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

######################### Make tables #########################

tables = []
for key in operations:
	op = operations[key]

	for regime in regimes:
		header = [regime+key] 
		header.extend(names)

		CFtable = []
		CFtable.append(Cnum)

		signal = cutflows["diHiggs/"+regime]
		for folder in folders:
			if folder=='diHiggs/' and ( key != '' and key != ' Nevt' ) :
				CFtable.append(len(Cnum)*['-'])
			else:
				row = op(signal,cutflows[folder+regime])
				CFtable.append(row)
                                
		CFtable=zip(*CFtable)
		tables.append(tabulate(CFtable, header,numalign="center", floatfmt=".1E", tablefmt=tablefmt))

######################### Print tables #########################

i=1
for table in tables:
        f = open('table_'+labeltag+"_"+str(i)+'.tex', 'w')
        i = i + 1
        for row in table:
                s=str(row)
                f.write(s)

        f.close()
