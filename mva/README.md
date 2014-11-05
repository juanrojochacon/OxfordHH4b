MVA README
==========

In this folder we have the code for the MVA analysis

This code takes the output of the HH4b code (in NTuples)
and trains an ANN via a naive genetic algorithm.

It requires that you have installed:

* GSL 

************************************

Input data
---------------------

MVA reads plaintext Ntuples produced by HH4b.
These may be for any kinematic variable choices.

The code takes the nTuples from it's first argument e.g

./mva ../analysis/res/ucl/total/ntuple.dat

************************************

Output and plotting
---------------------

The code outputs the results of the MVA classification
to ./res/<MVA>.dat where <MVA> details the NN architecture
and number of generations used in the fit.

ROC curves and signal/background identification histograms
may be obtained by running the script ./res/plotMVAdat.py
on these datafiles. ROC curves comparing two MVA results
can be obtained by passing multiple arguments to plotMVAdat.py.

### Network plotting

If you have the NetworkX package installed, you can plot the neural network output directly.
The plotNetwork.py script takes a .net file (output from mva) and plots the relative weights
of the inputs.

************************************

