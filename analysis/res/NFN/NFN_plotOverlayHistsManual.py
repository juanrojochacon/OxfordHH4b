#! /usr/bin/python
from matplotlib import pyplot as plt
from operator import add, sub
import sys, os, math
import numpy as np
from FIKRI_helpers import *

colours = ['r', 'b', 'g', 'm', 'c', 'y', 'k']

stringAna = "boost"


# toBePlotted = [
#   {
#     "title"     : stringAna + " Signal Region",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "2Tag diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
#       "4Tag diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
#     },
#     "pdfname"    : "Compare_2tag4tag_diHiggs_"+ stringAna +"_SR_mHH"
#   },
#   {
#     "title"     : stringAna + " Sideband Region",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "2Tag diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SDBA_2tag.dat",
#       "4Tag diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SDBA_4tag.dat",
#     },
#     "pdfname"    : "Compare_2tag4tag_diHiggs_"+ stringAna +"_SBA_mHH"
#   },
#   {
#     "title"     : stringAna + " 2 Tag",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "Sideband diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SDBA_2tag.dat",
#       "SR diHiggs"       : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
#     },
#     "pdfname"    : "Compare_SidebandSig_diHiggs_"+ stringAna +"_2tag_mHH"
#   },
#   {
#     "title"     : stringAna + " 4 Tag",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "Sideband diHiggs" : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SDBA_4tag.dat",
#       "SR diHiggs"       : "baseline_noPU_atlas_qcd/diHiggs/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
#     },
#     "pdfname"    : "Compare_SidebandSig_diHiggs_"+ stringAna +"_4tag_mHH"
#   },
#   {
#     "title"     : stringAna + " Signal Region",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "2Tag Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
#       "4Tag Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
#     },
#     "pdfname"    : "Compare_2tag4tag_background_"+ stringAna +"_SR_mHH"
#   },
#   {
#     "title"     : stringAna + " Sideband Region",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "2Tag Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SDBA_2tag.dat",
#       "4Tag Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SDBA_4tag.dat",
#     },
#     "pdfname"    : "Compare_2tag4tag_background_"+ stringAna +"_SBA_mHH"
#   },
#   {
#     "title"     : stringAna + " 2 Tag",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "Sideband Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SDBA_2tag.dat",
#       "SR Bkgd"       : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SIG_2tag.dat",
#     },
#     "pdfname"    : "Compare_SidebandSig_background_"+ stringAna +"_2tag_mHH"
#   },
#   {
#     "title"     : stringAna + " 4 Tag",
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "Sideband Bkgd" : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SDBA_4tag.dat",
#       "SR Bkgd"       : "baseline_noPU_atlas_qcd/background/histo_m_HH_zoom_"+ stringAna +"_SIG_4tag.dat",
#     },
#     "pdfname"    : "Compare_SidebandSig_background_"+ stringAna +"_4tag_mHH"
#   },
# ]

# doNormUnity = False
# doLogY = False
# doLogX = False


AnaGoodName    = {"res"  : "Resolved", "boost" : "Boosted"}
TagGoodName    = {"2tag" : "2-Tag", "4tag" : "4-Tag"}
RegionGoodName = {"SIG"  : "Signal Region", "SDBA" : "Sideband"}
PlotGoodName   = {
"m_HH"  : "Mass(HH) [GeV]", 
"m_H0"  : "Mass(Lead H) [GeV]", 
"m_H1"  : "Mass(SubLead H) [GeV]",
"pt_H0" : "pT(Lead H) [GeV]", 
"pt_H1" : "pT(SubLead H) [GeV]",
"eta_H0" : "Eta(Lead H)", 
"eta_H1" : "Eta(SubLead H)"
}

stringAna    = "res"
stringRegion = "SIG" 
stringPlot   = "eta_H1"

toBePlotted = [
    {
    "title"     : AnaGoodName[stringAna] + " " + RegionGoodName[stringRegion],
    "xaxisname" : PlotGoodName[stringPlot] ,
    "plots"     :
    {
      "2-Tag Bkgd"   : "baseline_noPU_atlas_qcd/background/histo_"+stringPlot+"_"+ stringAna +"_"+stringRegion+"_2tag.dat",
      "4-Tag Bkgd"   : "baseline_noPU_atlas_qcd/background/histo_"+stringPlot+"_"+ stringAna +"_"+stringRegion+"_4tag.dat",
    },
    "pdfname"    : "Compare_2tag4tag_bkgd_"+ stringAna +"_"+stringRegion+"_"+stringPlot
  },
]
doNormUnity = True
doLogY = False
doLogX = False

# stringAna    = "res"
# stringTag    = "4tag" 
# stringRegion = "SIG" 
# toBePlotted = [
#   {
#     "title"     : AnaGoodName[stringAna] + " " + TagGoodName[stringTag] + " " + RegionGoodName[stringRegion],
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "QCD 4b"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4b/histo_m_HH_zoom_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#     },
#     "pdfname"    : "Histo_" + stringTag + "_QCD4b_"+ stringAna +"_"+stringRegion+"_mHH"
#   },
#   {
#     "title"     : AnaGoodName[stringAna] + " " + TagGoodName[stringTag] + " " + RegionGoodName[stringRegion],
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "QCD 2b2j"   : "baseline_noPU_atlas_qcd/SHERPA_QCD2b2j/histo_m_HH_zoom_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#     },
#     "pdfname"    : "Histo_" + stringTag + "_QCD2b2j_"+ stringAna +"_"+stringRegion+"_mHH"
#   },
#   {
#     "title"     : AnaGoodName[stringAna] + " " + TagGoodName[stringTag] + " " + RegionGoodName[stringRegion],
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "QCD 4j"   : "baseline_noPU_atlas_qcd/SHERPA_QCD4j/histo_m_HH_zoom_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#     },
#     "pdfname"    : "Histo_" + stringTag + "_QCD4j_"+ stringAna +"_"+stringRegion+"_mHH"
#   },
#   {
#     "title"     : AnaGoodName[stringAna] + " " + TagGoodName[stringTag] + " " + RegionGoodName[stringRegion],
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "ttbar"   : "baseline_noPU_atlas_qcd/SHERPA_QCDttbar/histo_m_HH_zoom_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#     },
#     "pdfname"    : "Histo_" + stringTag + "_QCDttbar_"+ stringAna +"_"+stringRegion+"_mHH"
#   },
# ]

# stringAna    = "res"
# stringTag    = "4tag" 
# stringRegion = "SDBA" 
# toBePlotted = [
#   {
#     "title"     : AnaGoodName[stringAna] + " " + TagGoodName[stringTag] + " " + RegionGoodName[stringRegion],
#     "xaxisname" : "Mass(HH) [GeV]",
#     "plots"     :
#     {
#       "QCD 4b"  : "baseline_noPU_atlas_qcd/SHERPA_QCD4b/histo_m_HH_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#       "QCD 2b2j": "baseline_noPU_atlas_qcd/SHERPA_QCD2b2j/histo_m_HH_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#       "QCD 4j"  : "baseline_noPU_atlas_qcd/SHERPA_QCD4j/histo_m_HH_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#       "ttbar"   : "baseline_noPU_atlas_qcd/SHERPA_QCDttbar/histo_m_HH_"+ stringAna +"_"+stringRegion+"_" + stringTag + ".dat",
#     },
#     "pdfname"    : "Overlay_" + stringTag + "_allbkgd_"+ stringAna +"_"+stringRegion+"_mHH"
#   },
# ]
# doNormUnity = False
# doLogY = False
# doLogX = False

for toPlot in toBePlotted:
  
  title     = toPlot["title"]
  xaxisname = toPlot["xaxisname"]
  plots     = toPlot["plots"]
  pdfname   = toPlot["pdfname"]

  # Setup figure
  fig, ax = plt.subplots()
  
  if(doLogY):
    ax.set_yscale('log',nonposx='clip')
  if(doLogX):
    ax.set_xscale('log',nonposy='clip')
    
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
      
    if doNormUnity:
      norm = np.sum(norm) # Numpy types: Normalize to same area 
    else:
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
    
  #plt.ylim(ymin=0.01)

  # Gridlines
  ax.xaxis.grid(True)
  ax.yaxis.grid(True)

  # Legend
  legend = ax.legend(loc='best')
  legend.get_frame().set_alpha(0.8)

  fig.savefig(pdfname + '.pdf')
