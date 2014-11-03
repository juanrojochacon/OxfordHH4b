#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy as np

colours = ['r', 'b', 'g']
icol = 0

# Setup figure
fig, ax = plt.subplots()
#ax.set_yscale('log')

ax.set_ylabel("Arbitary units")
ax.set_xlabel("Observable bin")

for idat in xrange(1,len(sys.argv)):
  
  infilenm = sys.argv[idat]
  basename = os.path.splitext(infilenm)[0]
    
  # Verify paths
  if os.path.exists(infilenm) == False:
    print "Error: source file" + infilenm + " not found!"
    sys.exit()
  
  infile = open(infilenm, 'rb')
  print "Processing " + infilenm + " ..."
  datafile = open(infilenm, 'rb')

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

  # Normalisations
  norm = 0
  lastbin = 0
  for i in xrange(0,len(xbin)):
    h = xbin[i] - lastbin
    lastbin = xbin[i]
    norm = norm + h*yval[i]

  norm = np.sum(norm) # Numpy types


  xlow = 0
  for x in xrange(0,len(xbin)):
    xvals = [xlow, xbin[x]]
    yvalsup = [CVup[x], CVup[x]]
    yvalsdn = [CVdn[x], CVdn[x]]
    xlow = xbin[x]
    ax.fill_between(xvals, yvalsup/norm, yvalsdn/norm, facecolor=colours[icol], alpha = 0.4, linewidth = 1, color = colours[icol])

  ax.plot(xbin,yval/norm,ls = "steps-pre", color = colours[icol], label=infilenm)
  icol=icol+1

# Gridlines
ax.xaxis.grid(True)
ax.yaxis.grid(True)

# Legend
legend = ax.legend(loc='best')
legend.get_frame().set_alpha(0.8)

fig.savefig('histo.pdf')
