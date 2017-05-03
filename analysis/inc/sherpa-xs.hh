#ifndef SHERPA_XS_H
#define SHERPA_XS_H
#include <string>
#include <vector>

int get_evt_info(std::string& tmp, double& wgt, double& wgt_norm, double& n_trials);
int get_xsec_info(std::string& tmp, double& gen_xsec);

std::vector<std::string>& split(const std::string& s, char delim, std::vector<std::string>& elems);

std::vector<std::string> split(const std::string& s, char delim);

#endif
