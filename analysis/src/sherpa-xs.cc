// =================================================================================== //
// Functions for calculating the cross-section from a HepMC2 file.
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// Example usage: (should just make standalone; could be a useful utility)
/* -----------------------------------------------------------------------------------
  double sum_wgts = 0., sum_trls = 0.; tot_xsec = 0.;\
  double wtt = 0., wnt = 0., ntt = 0.;
  int tmp = 0;

  ifstream hepmc_is ( ... );      // HepMC input

  // Loop through once to get the total trials and final cross-section
  while( hepmc_is ){
    getline( hepmc_is, line );
    tmp = get_evt_info( line, wtt, wnt, ntt );
    if( tmp ){
      sum_wgts += wtt;
      sum_trls += ntt;
    }
  }

  tot_xsec = sum_wgts/sum_trls;
-------------------------------------------------------------------------------------- */
//
//   @author Fady Bishara
//   @version 1.0 -- 01/04/2016
// =================================================================================== //

#include "sherpa-xs.hh"
#include <sstream>
#include <string>
#include <vector>

// Takes a line from a HepMC2 event file.
// -- Returns 1 if line starts with 'E'; i.e. marks the beginning of an event and so could
//    be used to count the number of events
// -- Returns 0 otherwise
int get_evt_info(std::string& tmp, double& wgt, double& wgt_norm, double& n_trials) {

    std::string evtstr = "E ";

    if (tmp.length() == 0) return 0;
    if (tmp.compare(0, evtstr.length(), evtstr) == 0) {
        std::vector<std::string> temp = split(tmp, ' ');
        if (temp.size() >= 17) {
            wgt      = std::stod(temp[13]);
            wgt_norm = std::stod(temp[13]) / std::stod(temp[15]);
            n_trials = std::stod(temp[16]);
        }
        return 1;
    }
    else
        return 0;
}

// Functions for splitting a string by a delimiter.
// Taken from: http://stackoverflow.com/a/236803
// with the improvement suggested in one of the comments to push only non-empty elements
std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems) {
    std::stringstream ss(s);
    std::string       item;
    while (std::getline(ss, item, delim)) {
        if (!item.empty()) elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string& s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
