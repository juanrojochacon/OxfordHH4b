
// C++
#include <fstream>
#include <unistd.h>
#include <dirent.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenCrossSection.h"
#include "HepMC/IO_GenEvent.h"

// Merge HepMC samples
// This merging procedure is /VERY SHAKY/ for unweighted events.

std::vector<std::string> getHepMC( std::string pth )
{
    std::vector<std::string> files;
    const char* workingDirectory = pth.c_str();
    const auto dir = opendir(workingDirectory);
    if (dir) {
        auto dirent = readdir(dir);
        while (dirent) {
            if (dirent->d_type != DT_DIR) {
                std::string name{dirent->d_name, dirent->d_namlen};
                if (name.find(".hepmc") != std::string::npos) 
                  files.push_back( pth +'/'+ name );
            }
            dirent = readdir(dir);
        }
        closedir(dir);
    }
    return files;
}

int main( int argc, char* argv[] ) 
{  
 if (argc != 2)
  {
    std::cerr << "Error: Wrong number of arguments!"<<std::endl;
    std::cerr << "Usage: mergeHepMC <directory for samples>" <<std::endl;
    std::cerr << "This code will merge all .hepmc files in the target directory" <<std::endl;
    exit(-1);
  }

  int counter = 0; double sum_weights = 0; double sum_trials = 0;
  HepMC::IO_GenEvent ascii_out("merge.hepmc",std::ios::out);

  const auto files = getHepMC(argv[1]);
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
