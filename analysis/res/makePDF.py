#!/usr/bin/python

from pdfMakerFromPDF import PDFMakerFromEPS
import os

#print argv

def mySnazzyPdfMaker(test):
    print test


fileType = '.pdf'

plotPath = os.path.join("/home","frostj","OxfordHbb","oxfordcutsSVN","OxfordHH4b","branches","oxfordcuts","analysis","res","boost_normplots")
paths = []
for file in os.listdir(plotPath):
    if file.endswith(fileType):
        paths.append(os.path.join(plotPath, file))

Maker = PDFMakerFromEPS(paths, 'boost_norm.pdf', fileType)
Maker.writePDF()

