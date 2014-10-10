**************
README
*************

In this folder we have the main analysis code 

HH4b.cc

which takes a Les Houches parton level event file, showers it with
Pythia8 and then does the corresponding analysis, filling histograms,
creating plots etc.

It requires that you have installed:

* Pythia8 (I have tested up to v8.186):

* Root

* Fastjet
http://fastjet.fr/

Here one also needs to install fastjet contrib
http://fastjet.hepforge.org/contrib/
to access the Variable R jets


* LHAPDFv5

The config programs for each of these dependencies should be reachable
by the makefile. You should edit the Makefile samples path to the
Dropbox location

************************************

To compile the variableR example code, type

make variableR

************************************

