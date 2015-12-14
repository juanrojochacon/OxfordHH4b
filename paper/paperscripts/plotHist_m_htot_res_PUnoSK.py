#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy as np

datafiles = []
plotnames = []

############################## Settings ###############################

# Output filename
outfile = "m_htot_res_signal_PUnoSK"

linestyles = [ 'solid', 'dashed','dotted', 'dashed']

# Datafiles
datafiles.append("../plotdata/20K_noPU/histo_m_H_res_C2.0.dat")
datafiles.append("../plotdata/20K_PU80SK/histo_m_H_res_C2.0.dat")
datafiles.append("../plotdata/20K_PU80/histo_m_H_res_C2.0.dat")


# Plot labels
plotnames.append("no PU")
plotnames.append(r"$\langle n_{PU} \rangle =80$ + SK")
plotnames.append(r"$\langle n_{PU} \rangle =80$")


# Axis labels
xLabel = "Higgs Candidates Invariant Mass (GeV)"
yLabel = "a. u."

# Log axes
xLog = False
yLog = False

# Normalise histograms
normalise = True

#######################################################################

if len(datafiles)!=len(plotnames):
  print "Error: datafile and plotname arrays are different lengths!"
  exit()

colours = ['r', 'b', 'k']
icol = 0

# Setup figure
fig, ax = plt.subplots()
if xLog == True:
  ax.set_xscale('log')
if yLog == True:
  ax.set_yscale('log')

ax.set_ylabel(yLabel)
ax.set_xlabel(xLabel)


ax.set_xlim([85,170])
ax.set_ylim([0,0.06])




for idat in xrange(0,len(datafiles)):

  # Verify paths
  infilenm = datafiles[idat]
  if os.path.exists(infilenm) == False:
    print "Error: source file" + infilenm + " not found!"
    sys.exit()
  
  infile = open(infilenm, 'rb')
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
  norm = 1.0
  if normalise == True:
    norm = 0
    for i in xrange(0,len(xhi)):
      h = xhi[i] - xlo[i]
      norm = norm + h*yval[i]
  norm = np.sum(norm) # Numpy types

  
  # Error bars
  CVup = map(add, yval/norm, errup/norm)
  CVdn = map(sub, yval/norm, errdn/norm)

  for x in xrange(0,len(xhi)):
    xvals = [xlo[x], xhi[x]]
    yvalsup = [CVup[x], CVup[x]]
    yvalsdn = [CVdn[x], CVdn[x]]

  # Insert lower x-values
  xhi.insert(0,xlo[0])
  yval.insert(0,yval[0])

  ax.plot(xhi,yval/norm,drawstyle = "steps-pre", color = colours[icol], label=plotnames[idat],linestyle=linestyles[idat],linewidth=2.4)
  icol=icol+1

plt.rcParams.update({'font.size': 14})
fig.text(0.21,0.93,"Signal events, Resolved category", fontsize=20)

# Legend
legend = ax.legend(loc='best')
legend.get_frame().set_alpha(0.8)

fig.savefig(outfile+'.pdf')
