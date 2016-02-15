//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

// Modified for OxfordHH4b - nh 22/10/14

#pragma once

#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

/**
 *  \class Parametrization
 *  \brief Virtual parametrizaton base class
 */
class Parametrization
{
public:
  Parametrization(string name);
  Parametrization(Parametrization const& );
  virtual ~Parametrization();
  
  virtual void Compute(const double*,double*) const = 0;

  void CopyPars(Parametrization*);
  virtual Parametrization* Duplicate() = 0;
  
  double* GetParameters() {return fParameters;}
  int const& GetNParameters() const {return fNParameters;}
  
  string const& GetParamName() const {return fParamName;};
  
  void ExportPars(string const& path) const;
  void ImportPars(string const& path);

protected:
  int   fNParameters;   //!< Total number of parameters
  double* fParameters;    //!< Parameters
  
private:
  const string fParamName;

};

/**
 *  \class MultiLayerPerceptron
 *  \brief General MLP class
 */

class MultiLayerPerceptron : public Parametrization
{
public:
  MultiLayerPerceptron(vector<int> const& arch);         //!< Network constructor
  MultiLayerPerceptron(MultiLayerPerceptron const&);  //!< Network copy constructor
  ~MultiLayerPerceptron();                            //!< Network destructor
  
  void Compute(const double*,double*) const;  //!< Returns a fArch[fNLayers-1] long array of output for a given input array
  Parametrization* Duplicate();     //!< Returns a parametrization based on an MLP
  
  const int*  GetArch() const {return fArch;};
  const int   GetNumNodeParams(int const& layer);                   //!< Returns the number of parameters that belong to a specific node (including biases).
  double*       GetNodeParams   (int const& layer, int const& node);  //!< Returns a pointer to the fParameters coordinate representing the parameters for a specific node

  void Export(string const& path, string const& kinstream) const;
    
private:
  const int fNLayers;   //!< Number of layers
  int* fArch;           //!< Network architecture
  
  double** fWeightMatrix;         //!< Weights/Thresholds
  mutable double** fOutputMatrix; //!< Neuron Activation
};

