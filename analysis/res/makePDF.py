#!/usr/bin/python

from pdfMakerFromPDF import PDFMakerFromEPS
import os

#print argv

def mySnazzyPdfMaker(test):
    print test


fileType = '.pdf'

plotPath = os.path.join("/home","behr","Private","OxfordHH4b.git","trunk","analysis","res","plots")
paths = []
for file in os.listdir(plotPath):
    if file.endswith(fileType):
        paths.append(os.path.join(plotPath, file))

Maker = PDFMakerFromEPS(paths, 'boost_unnorm.pdf', fileType)
Maker.writePDF()

