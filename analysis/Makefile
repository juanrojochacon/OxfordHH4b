SHELL = /bin/sh

###################################
####################################
####################################
#
# USER DEFINED PATHS
#
# Path to Pythia8 installation
PYTHIA8DIR=/Users/juanrojo/physics12/HHvbf/pythia8180
#
# Path to LHAPDF libraries
LHAPDFLOCATION=/Users/juanrojo/software/lhapdf/lha-5.9.0/lib
#
# Path to ROOT installation
ROOTDIR=/Users/juanrojo/software/root/root-v5-34-00-patches
#
# Path to FastJet installation
FASTJETPATH=/Users/juanrojo/software/fastjet/fastjet-3.0.3-install
#
# The rest of the Makefile should not be edited
#
####################################
####################################
####################################

# Makefile options
-include $(PYTHIA8DIR)/config.mk

# Location of directories.
INCDIR=$(PYTHIA8DIR)/include
LIBDIRARCH=$(PYTHIA8DIR)/lib/archive
BINDIR=bin

# root libraries
ROOTLIBS=-L$(ROOTDIR)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,$(ROOTDIR)/lib -lm -ldl

# fastjet libraries
FASTJETLIBS=-Wl,-rpath,$(FASTJETPATH)/lib -lm -L$(FASTJETPATH)/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone

# Create an executable that links to LHAPDF
HH4b variableR: $(PYTHIA8LOCATION)/$(LIBDIRARCH)/libpythia8.a
	@mkdir -p $(BINDIR)
	$(CXX) $(CXXFLAGS) -I$(ROOTDIR)/include -I$(PYTHIA8LOCATION)/$(INCDIR) $@.cc \
	-I$(FASTJETPATH)/include -o $(BINDIR)/$@.exe -L$(PYTHIA8LOCATION)/$(LIBDIRARCH) -lpythia8  \
	-L$(LHAPDFLOCATION) -lLHAPDF $(ROOTLIBS) \
	-L$(FASTJETPATH)/lib -lVariableR -L$(FASTJETPATH)/lib $(FASTJETLIBS) $(FLIBS)
	@ln -fs $(BINDIR)/$@.exe $@.exe

# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf $(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core* *.eps *.root *.C *.dat *.res *.tex *.pdf

#####################################################################################################
