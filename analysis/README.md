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

************************************

Plotting results
******************
