#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math

colours = ['y', 'r', 'g']
icol = 1

if len(sys.argv) < 2:
	print "Error: please specify a histogram file!"
	sys.exit()

# Plot source
source = sys.argv[1]

# Verify path
if os.path.exists(source) == False:
  	print "Error: source file" + source + " not found!"
  	sys.exit()

# Open commondata
print "Plotting: " + source
datafile = open(source, 'rb')

xbin = []
yval = []

errup = []
errdn = []

dataread = False
for line in datafile:
	linesplit = line.split()

	if len(linesplit) == 1:
		continue

	if linesplit[1] == 'END':
		break

	if dataread == True:
		xbin.append(float(linesplit[1]))
		#xbin.append(float(linesplit[0]))
		yval.append(float(linesplit[2]))
		errdn.append(float(linesplit[3]))
		errup.append(float(linesplit[4]))

	if linesplit[1] == 'xlow':
		dataread = True

# Error bars
CVup = map(add, yval, errup)
CVdn = map(sub, yval, errdn)

fig, ax = plt.subplots()
#if min(yval) == 0:
#	if math.log(max(yval)) > 100:
#		ax.set_yscale('log')
#elif max(yval)/min(yval) > 100:
#	ax.set_yscale('log')
ax.set_yscale('log')

xlow = 0
for x in xrange(0,len(xbin)):
	xvals = [xlow, xbin[x]]
	yvalsup = [CVup[x], CVup[x]]
	yvalsdn = [CVdn[x], CVdn[x]]
	xlow = xbin[x]
	ax.fill_between(xvals, yvalsup, yvalsdn, facecolor=colours[icol], alpha = 0.4, linewidth = 1, color = colours[icol])

ax.plot(xbin,yval,ls = "steps-pre", color = colours[icol])
fig.savefig(os.path.splitext(source)[0]+'.pdf')
fig.savefig(os.path.splitext(source)[0]+'.png', dpi=80)
