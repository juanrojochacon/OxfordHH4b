#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy as np
from FIKRI_helpers import *

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']

stringAna = "res"


toBePlotted = [
  {
    "title"     : stringAna + " 4-tag Signal Region",
    "xaxisname" : "Mass(HH) [GeV]",
    "plots"     :
    {
      "QCD 4b"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4b/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
      "QCD 2b2j" : "baseline_noPU_atlas_qcd/SHERPA_QCD2b2j/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
      "QCD 4j"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4j/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
      "ttbar"    : "baseline_noPU_atlas_qcd/SHERPA_QCDttbar/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
    },
    "pdfname"    : "Stack_4tag_backgrounds_"+ stringAna +"_SR_mHH"
  },
  {
    "title"     : stringAna + " 2-tag Signal Region",
    "xaxisname" : "Mass(HH) [GeV]",
    "plots"     :
    {
      "QCD 4b"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4b/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
      "QCD 2b2j" : "baseline_noPU_atlas_qcd/SHERPA_QCD2b2j/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
      "QCD 4j"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4j/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
      "ttbar"    : "baseline_noPU_atlas_qcd/SHERPA_QCDttbar/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
    },
    "pdfname"    : "Stack_2tag_backgrounds_"+ stringAna +"_SR_mHH"
  },
]

for toPlot in toBePlotted:
  
  title     = toPlot["title"]
  xaxisname = toPlot["xaxisname"]
  plots     = toPlot["plots"]
  pdfname   = toPlot["pdfname"]

  # Setup figure
  fig, ax = plt.subplots()
  #ax.set_yscale('log')
  ax.set_xlabel(xaxisname)
  ax.set_ylabel("Arbitary units")
  plt.title(title)
  
  icol = 0
  for plot in plots:
    
    infilenm  = plots[plot] 

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

    #norm = np.sum(norm) # Numpy types: Normalize to same area 
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
    xhi.insert(0, xlo[0])
    yval.insert(0,yval[0])

    ax.plot(xhi, yval/norm, ls = "steps-pre", color = colours[icol], label=plot)
    icol=icol+1

  # Gridlines
  ax.xaxis.grid(True)
  ax.yaxis.grid(True)

  # Legend
  legend = ax.legend(loc='best')
  legend.get_frame().set_alpha(0.8)

  fig.savefig(pdfname + '.pdf')
