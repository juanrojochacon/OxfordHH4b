#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math, numpy
from pylab import pcolor, show, colorbar, xticks, yticks

if len(sys.argv) < 2:
	print "Error: please specify an NTuple file!"
	sys.exit()

# Ntuple source
ntup = sys.argv[1]
datafile = open(ntup, 'rb')

line = datafile.readline()
nKin = len(line.split()) - 4
print nKin, "Variables Read"

# Read kinematics labels
kinLabels = []
for i in xrange(4,nKin+4):
	kinLabels.append(line.split()[i])

# Read kinematics
kinematics = []
for i in xrange(0,nKin):
	kinematics.append([])

for line in datafile:
	linesplit = line.split()

	for i in xrange(0,nKin):
		kinPoint = float(linesplit[i+3])
		kinematics[i].append(kinPoint)

# Calculate covariance matrix
covMat = numpy.corrcoef(kinematics) 

# Plot correlation matrix
plt.xticks(rotation=70)
fig, ax = plt.subplots()
cax = ax.pcolor(covMat, vmin=-1, vmax=1)

# Adjustments
ax.axis([0,9,0,9])
fig.subplots_adjust(bottom=0.15, left=0.15, right=1.02)

# Colorbar
fig.colorbar(cax)

# Title
fig.text(0.15,0.92,os.path.basename(ntup), fontsize=12)

# put the major ticks at the middle of each cell
ax.set_xticks(numpy.arange(covMat.shape[0])+0.5, minor=False)
ax.set_yticks(numpy.arange(covMat.shape[1])+0.5, minor=False)

ax.set_xticks(numpy.arange(covMat.shape[0]), minor=True)
ax.set_yticks(numpy.arange(covMat.shape[0]), minor=True)

# Set ticks
ax.set_xticklabels(kinLabels, minor=False, rotation=-30)
ax.set_yticklabels(kinLabels, minor=False)

# Gridlines
ax.xaxis.grid(True, which = 'minor')
ax.yaxis.grid(True, which = 'minor')

fig.savefig("correlation.pdf")
