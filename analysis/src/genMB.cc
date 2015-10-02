// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Mikhail Kirsanov, Mikhail.Kirsanov@cern.ch, based on main01.cc.
// This program illustrates how HepMC can be interfaced to Pythia8.
// It studies the charged multiplicity distribution at the LHC.
// HepMC events are output to the hepmcout41.dat file.

// WARNING: typically one needs 25 MB/100 events at the LHC.
// Therefore large event samples may be impractical.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

int main( int argc, char* argv[]  ) {


 if (argc != 2)
  {
    cerr << "Error: Wrong number of arguments!"<<endl;
    cerr << "Usage: HH4b <runID>" <<endl;
    exit(-1);
  }

  // Read sampleID
  const int runID = atoi(argv[1]);

  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC::Pythia8ToHepMC ToHepMC;

  // Specify file where HepMC events will be stored.
  std::stringstream filename;
  filename << "/data/atlastemp/DiHiggsSharedSamples/PYTHIA_MinBias_14TEV." << runID <<".hepmc";
  std::cout << "Using filename: "<<filename.str()<<endl;
  HepMC::IO_GenEvent ascii_io(filename.str(), std::ios::out);

  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia(std::string(PYTHIADIR));
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  pythia.readString("Random:setSeed = on");

  std::stringstream seedString;
  seedString << "Random:seed = " <<  ( 23487*runID + 4847 );
  pythia.readString(seedString.str());
  pythia.readString("SoftQCD:all = on");
  pythia.settings.parm("Beams:eCM", 14000);
  pythia.init();

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < 1E7; ++iEvent) {
    if (!pythia.next()) continue;

    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build; but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop. Statistics. Histogram.
  }
  pythia.stat();

  // Done.
  return 0;
}
