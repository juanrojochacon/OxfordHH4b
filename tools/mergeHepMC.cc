
// C++
#include <fstream>
#include <unistd.h>
#include <dirent.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/IO_GenEvent.h"

// Merge HepMC samples
// This merging procedure is /VERY SHAKY/ for unweighted events.

int main( int argc, char** argv ) 
{  
 if (argc < 2)
  {
    std::cerr << "Error: Wrong number of arguments!"<<std::endl;
    std::cerr << "Usage: mergeHepMC <samples>" <<std::endl;
    std::cerr << "This code will merge all .hepmc files in arguments" <<std::endl;
    exit(-1);
  }

  int counter = 0; double sum_weights = 0; double sum_trials = 0;
  HepMC::IO_GenEvent ascii_out("merge.hepmc",std::ios::out);

  std::vector<std::string> files(argv + 1, argv + argc);
  for (auto file : files )
  {
    std::cout << "Reading: " <<file <<std::endl;
    HepMC::IO_GenEvent ascii_in(file.c_str(),std::ios::in);
    HepMC::GenEvent* evt =  ascii_in.read_next_event();
    while ( evt ) {
      counter++;
      sum_weights += evt->weights()[0];
      sum_trials += evt->weights()[3];

      if (counter % 10000 == 0) 
      {
        const double gen_xsec = evt->cross_section()->cross_section();
        std::cout << "Processed " <<counter<<" events"<<std::endl;
        std::cout << gen_xsec <<"  "<<sum_weights/sum_trials<<std::endl;
      }

      evt->set_event_number(counter);
      evt->cross_section()->set_cross_section(sum_weights/sum_trials);

      ascii_out << evt;
      delete evt;
      ascii_in >> evt;
    }
  }

  // End of the main program
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
