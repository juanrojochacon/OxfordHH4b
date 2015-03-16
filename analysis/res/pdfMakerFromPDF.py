import os

class PDFMakerFromEPS(object):
  def __init__(self, epsPathsAndFrameTitles, saveName = 'plots.pdf', imageType = '.eps'):
    self.imageType = imageType
    self._epsPathToFrameTitles = self._generateEPSPathToFrameTitlePairs(epsPathsAndFrameTitles)

    self.pdfFilePath = saveName
    #self.imageType = imageType
    self.writeTitlePage = True
    self.title = 'Plot collection'
    self.author = 'K Behr'

    self._latexFileName = ''
    self._latexFile = None

  def writePDF(self):
    self._createLatexFile()
    self._compileLatexFile()
    self._renamePDF()
    self._deleteIntermediateFiles()

  def _generateEPSPathToFrameTitlePairs(self, epsPathsAndFrameTitles):
    import os.path
    if type(epsPathsAndFrameTitles[0]) == str:
      pathToTitlePairs = []
      for path in epsPathsAndFrameTitles:
        frameTitle = os.path.basename(path)
        frameTitle = frameTitle.replace(self.imageType, '')
        frameTitle = frameTitle.replace('_', '-')
        pathToTitlePairs.append((path, frameTitle))
      return pathToTitlePairs
    elif type(epsPathsAndFrameTitles[0]) == tuple:
      return epsPathsAndFrameTitles
  
  def _getWorkingFolder(self):
    import os.path
    pdfDirectory = os.path.dirname(self.pdfFilePath)
    if pdfDirectory == '':
      pdfDirectory = '.'
    return pdfDirectory

  def _absoluteDirname(self, path):
    absolutePath = os.path.abspath(path)
    dirname = os.path.dirname(absolutePath)
    return dirname

  def _pdfFileName(self):
    return os.path.basename(self.pdfFilePath)

  def _createLatexFile(self):
    self._latexFileName = self._getUnusedTemporaryName('tex')
    print '_latexFileName:',self._latexFileName
    self._latexFile = open(self._latexFileName, 'w')
    try:
      self._writeLatexHeader()
      self._writeOneFramePerPlot()
      self._writeLatexFooter()
    finally:
      self._latexFile.close()
      self._latexFile = None

  def _getUnusedTemporaryName(self, extension):
    fileName = self._pdfFileName().replace(self.imageType, '')
    while os.path.exists(os.path.join(self._getWorkingFolder(), fileName)+'.'+extension):
      fileName += '1'
    return os.path.join(self._getWorkingFolder(), fileName)+'.'+extension

  def _writeLatexHeader(self):
    self._latexFile.write('\\documentclass{beamer}\n')
    self._latexFile.write('\\usepackage{graphicx}\n')
    self._latexFile.write('\\usepackage{amssymb}\n')
    self._latexFile.write('\\usepackage{epstopdf}\n')
    self._latexFile.write('\\setbeamertemplate{footline}{\insertframenumber/\inserttotalframenumber}\n')
    self._latexFile.write('\\graphicspath{ {./} }\n')
    self._latexFile.write('\\title{'+self.title+'}\n')
    self._latexFile.write('\\author{'+self.author+'}\n')
    self._latexFile.write('\\institute{University of Oxford}\n')
    self._latexFile.write('\\begin{document}\n')
    if self.writeTitlePage:
      self._latexFile.write('\\frame{\\titlepage}\n')

  def _writeOneFramePerPlot(self):
    for epsFilePathAndTitle in self._epsPathToFrameTitles:
      epsFilePath = epsFilePathAndTitle[0]
      frameTitle = epsFilePathAndTitle[1]
      self._latexFile.write('\\frame{\n')
      self._latexFile.write('  \\frametitle{'+frameTitle+'}\n')
      self._latexFile.write('  \\begin{figure}\n')
      self._latexFile.write('    \\includegraphics[width=\\textheight,height=0.8\\textheight,keepaspectratio]{'+epsFilePath+'}\n')
      self._latexFile.write('  \\end{figure}\n')
      self._latexFile.write('}\n')

  def _writeLatexFooter(self):
    self._latexFile.write('\\end{document}\n')

  def _compileLatexFile(self):
    nameStem = self._latexFileNameStem()
    print 'Name stem:', nameStem
    print 'Working folder:', self._getWorkingFolder()
    command = 'cd '+self._getWorkingFolder()+';'
    if self.imageType == ".eps":
      command += 'latex {0};'.format(nameStem)
      command += 'latex {0};'.format(nameStem)
      command += 'dvips {0}.dvi -o {0}.ps;'.format(nameStem)
      command += 'ps2pdf {0}.ps {0}.pdf;'.format(nameStem)
      command += 'cd '+os.getcwd()
    elif self.imageType == ".pdf":
      command += 'pdflatex {0};'.format(nameStem);
      command += 'pdflatex {0};'.format(nameStem);
    else:
      raise("Unrecognised image type")

    import subprocess
    subprocess.call(command, shell=True)

  def _latexFileNameStem(self):
    nameStem = os.path.basename(self._latexFileName).replace('.tex', '')
    return nameStem

  def _renamePDF(self):
    originalFileName = self._latexFileName.replace('.tex', '')+'.pdf'
    print 'RenamePDF'
    print 'originalFileName:',originalFileName
    print 'pdfFilePath:',self.pdfFilePath
    import shutil
    shutil.move(originalFileName, self.pdfFilePath)

  def _deleteIntermediateFiles(self):
    import os
    junkExtensions = ['aux', 'dvi', 'log', 'nav', 'out', 'ps', 'snm', 'tex', 'toc']
    for extension in junkExtensions:
      junkFileName = os.path.join(self._getWorkingFolder(), self._latexFileNameStem())+'.'+extension
      if os.path.exists(junkFileName):
        os.remove(junkFileName)
