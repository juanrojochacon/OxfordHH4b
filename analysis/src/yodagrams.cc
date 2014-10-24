// yodagrams.cc 17/10/14 nh

#include "yodagrams.h"
//#include "settings.h"
#include "YODA/Histo1D.h"
#include "YODA/WriterFLAT.h"

// Static map of histograms
static std::map<int,YODA::Histo1D*> yodaHistograms;

  // General string hasher
static int IntHash(const std::string& _str)
{
	const char* s = _str.c_str();
	unsigned h = 31;
	while (*s) {
		h = (h * 54059) ^ (s[0] * 76963);
		s++;
	}
	return h % 86969;
};

// Basic hashing for stl maps 
class HashFunctor
{
	public: size_t operator()(const char* s) {
		return IntHash(s);
	}
};

void yoda_add(YODA::Histo1D* hist){

	std::map<int,YODA::Histo1D*>::iterator iMap = yodaHistograms.find(IntHash(hist->path()));
	if (iMap != yodaHistograms.end())
	{
		std::cerr << "yoda_add error: HASH COLLISION for histogram: "<<hist->path()<<std::endl;
		std::cerr << "Either histogram is duplicated, or we need a better hash function!"<<std::endl;
	}
	else
	{
		yodaHistograms.insert(std::make_pair(IntHash(hist->path()),hist));
	}

	hist->setAnnotation(std::string("XLabel"), std::string("p_T (GeV)"));
}

void yoda_fill(std::string const& histopath, double const& weight, double const& coord)
{
	std::map<int,YODA::Histo1D*>::iterator iMap = yodaHistograms.find(IntHash(histopath));
  if (iMap != yodaHistograms.end())
  {	(*iMap).second->fill(coord,weight);	}
  else
  {
    std::cerr << "yoda_fill error: Cannot find Histogram: "<<histopath<<std::endl;
    exit(-1);
  }
  
}

void yoda_export()
{
	std::map<int,YODA::Histo1D*>::iterator iMap = yodaHistograms.begin();
	while (iMap != yodaHistograms.end())
	{
		std::cout << "Writing Histogram: "<< (*iMap).second->path()<<std::endl;
		YODA::WriterFLAT::write("." + (*iMap).second->path(), *(*iMap).second);
		iMap++;
	}
}
