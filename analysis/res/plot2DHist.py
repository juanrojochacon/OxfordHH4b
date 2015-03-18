#! /usr/bin/python
from matplotlib import pyplot as plt
from collections import OrderedDict
import numpy as np
import sys, math, os

nArgs = len(sys.argv)
if nArgs != 2:
	print("plot2DHist only supports one argument!")
	sys.exit()

infilenm = sys.argv[1]
basename = os.path.splitext(infilenm)[0]

# Verify paths
if os.path.exists(infilenm) == False:
	print "Error: source file" + infilenm + " not found!"
	sys.exit()

infile = open(infilenm, 'rb')
datafile = open(infilenm, 'rb')

# Histogram bin edges
xEdges = []
yEdges = []

# Histogram data
xVals = []
yVals = []
weights = []

dataread = False
for line in datafile:
	linesplit = line.split()

	if len(linesplit) == 1:
	  continue

	if linesplit[1] == 'END':
	  break

	if dataread == True:
		xEdges.append(float(linesplit[0]))
		xEdges.append(float(linesplit[1]))
		yEdges.append(float(linesplit[2]))
		yEdges.append(float(linesplit[3]))

		xPos = (xEdges[-2] + xEdges[-1])/2.0
		yPos = (yEdges[-2] + yEdges[-1])/2.0

		weight = float(linesplit[4])

		xVals.append(xPos)
		yVals.append(yPos)
		weights.append(weight)
		
	if linesplit[1] == 'xlow':
	  dataread = True

# Remove duplicates
xEdges = list(OrderedDict.fromkeys(xEdges))
yEdges = list(OrderedDict.fromkeys(yEdges))

# Make numpy Histogram
H, xEdges, yEdges = np.histogram2d(xVals, yVals, bins=(xEdges, yEdges), normed=True, weights=weights)

# Flip (numpy returns the wrong ordering) and mask
H = np.rot90(H)
H = np.flipud(H)
Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
 
# Plot 2D histogram using pcolor
fig, ax = plt.subplots()
cax = ax.pcolormesh(xEdges,yEdges,Hmasked)

ax.xaxis.grid(True, which = 'majorminor')
ax.yaxis.grid(True, which = 'majorminor')

cbar = fig.colorbar(cax)
cbar.ax.set_ylabel('xSec [fb]')

fig.text(0.15,0.92,infilenm, fontsize=12)

fig.savefig('histo2D.pdf')