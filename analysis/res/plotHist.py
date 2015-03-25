#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy as np

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']
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

  xhi = []
  xlo = []
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
      xlo.append(float(linesplit[0]))
      xhi.append(float(linesplit[1]))
      yval.append(float(linesplit[2]))
      errdn.append(float(linesplit[3]))
      errup.append(float(linesplit[4]))

    if linesplit[1] == 'xlow':
      dataread = True

  # Normalisations
  norm = 0
  for i in xrange(0,len(xhi)):
    h = xhi[i] - xlo[i]
    norm = norm + h*yval[i]

  #norm = np.sum(norm) # Numpy types
  norm = np.sum(1.0) # Numpy types

  # Error bars
  CVup = map(add, yval/norm, errup/norm)
  CVdn = map(sub, yval/norm, errdn/norm)

  for x in xrange(0,len(xhi)):
    xvals = [xlo[x], xhi[x]]
    yvalsup = [CVup[x], CVup[x]]
    yvalsdn = [CVdn[x], CVdn[x]]
    ax.fill_between(xvals, yvalsup, yvalsdn, facecolor=colours[icol], alpha = 0.4, linewidth = 1, color = colours[icol])

  # Insert lower x-values
  xhi.insert(0,xlo[0])
  yval.insert(0,yval[0])

  ax.plot(xhi,yval/norm,ls = "steps-pre", color = colours[icol], label=infilenm[2:-4])
  icol=icol+1

# Gridlines
ax.xaxis.grid(True)
ax.yaxis.grid(True)

# Legend
legend = ax.legend(loc='best')
legend.get_frame().set_alpha(0.8)

fig.savefig('histo.pdf')
