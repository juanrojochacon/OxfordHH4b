SHELL = /bin/sh

###################################
####################################
####################################
#
# USER DEFINED PATHS

# Binary directory
BINDIR=bin
# MC sample directory
SAMPLEDIR= -DSAMPLEDIR=\"/Users/nathanhartland/Dropbox/HH4bMC/\"
# Pythia 8 xmldoc directory
PYTHIADIR = -DPYTHIADIR=\"$(shell pythia8-config --prefix)\"


# root libraries
ROOTFLAGS = $(shell root-config --cflags) $(shell root-config --glibs)
# fastjet libraries
FASTJETFLAGS = $(shell fastjet-config --libs) $(shell fastjet-config --cxxflags)
# Pythia 8 
PYTHIAFLAGS = $(shell pythia8-config --cxxflags) $(shell pythia8-config --libs)
# LHAPDF
LHAPDFFLAGS = $(shell lhapdf-config --cppflags) $(shell lhapdf-config --ldflags)


# All flags
CXXFLAGS = $(ROOTFLAGS) $(FASTJETFLAGS) $(PYTHIAFLAGS) $(LHAPDFFLAGS) $(PYTHIADIR) $(SAMPLEDIR)

#
# The rest of the Makefile should not be edited
#
####################################
####################################
####################################

# Create an executable that links to LHAPDF
HH4b variableR:
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) $@.cc -o $(BINDIR)/$@.exe 
	@ln -fs $(BINDIR)/$@.exe $@.exe

# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf $(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core* *.eps *.root *.C *.dat *.res *.tex *.pdf

#####################################################################################################
