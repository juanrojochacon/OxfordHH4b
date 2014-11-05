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

* YODA 1.3.0 (https://yoda.hepforge.org)

* Fastjet http://fastjet.fr/

Here one also needs to install fastjet contrib
http://fastjet.hepforge.org/contrib/
to access the Variable R jets

The config programs for each of these dependencies should be reachable
by the makefile. You should edit the Makefile samples path to the location of your samples.

************************************

To compile and run the code:

make
./HH4b

This runs the analysis code for all the analyses specified in ./src/HH4b.cc.
The code will output all results to /res/<analysis>/<sample>
where <analysis> specified the name of the type of analysis, and <sample> refers to the source .lhe file.
Additionally the "total" "signal" and "background" samples are generated.

Results are presented as YODA FLAT histograms for the histograms specified in the analysis.
Additionally the kinematics for input to the MVA are outputed to an "ntuple.dat" file, on
a sample by sample basis.

************************************

Writing analyses

An example barebones analysis class can be found in /src/basic.cc and /inc/basic.h.
This may be expanded upon to make a new analysis. For use in the HH4b main code,
the new analysis should be added as a module in the Makefile, and added to the
list of analyses in HH4b.cc

************************************

Plotting results

To plot the histograms output to the /res/ directory, use the plotHisto.py script
located in /res.

It's arguments are the histograms to be compared. E.g
./plotHisto.py ./ucl/signal/histo_m4b.dat ./ucl/background/histo_m4b.dat

To plot a correlation matrix upon an nTuple output, use the correlationPlot.py script
located in /res.

It's arguments are an nTuple, e.g
./correlationPlot.py ./ucl/signal/ntuple.dat 

************************************
