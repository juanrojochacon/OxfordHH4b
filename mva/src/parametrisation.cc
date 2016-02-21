// $Id: parametrization.cc 1506 2014-01-15 11:26:09Z s0673800 $
//
// NNPDF++ 2012
//
// Authors: Nathan Hartland,  n.p.hartland@ed.ac.uk
//          Stefano Carrazza, stefano.carrazza@mi.infn.it
//          Luigi Del Debbio, luigi.del.debbio@ed.ac.uk

// Modified for OxfordHH4b - nh 22/10/14

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include "parametrisation.h"
#include "random.h"

// ******************** Base ******************************
/**
 * @brief Parametrization constructor
 */
 Parametrization::Parametrization( string name ):
 fNParameters(0),
 fParameters(0),
 fParamName(name)
 {
 }

/**
 * @brief Parametrization constructor
 * @param o the copy contructor
 */
 Parametrization::Parametrization( Parametrization const& o):
 fNParameters(o.fNParameters),
 fParameters(new double[fNParameters]),
 fParamName(o.fParamName)
 {
  for (int i=0; i<fNParameters; i++)
    fParameters[i] = o.fParameters[i];
  
  return;
}

/**
 * @brief Parametrization destructor
 */
 Parametrization::~Parametrization()
 {
 }

/**
 * @brief Parametrization::CopyPars
 * @param t
 */
 void Parametrization::CopyPars(Parametrization* t)
 {
  if (fNParameters != t->GetNParameters())
  {
    cerr << "Parametrization::CopyPars Error: number of parameters does not match: "<< fNParameters<<" vs "<<t->GetNParameters()<<endl;
    exit(-1);
  }
  
  for (int i=0; i<fNParameters; i++)
    fParameters[i] = t->fParameters[i];
  
  return;
}

  void Parametrization::ExportPars(string const& path) const
  {
    std::ofstream out(path);
    for (int i=0; i<fNParameters; i++)
      out << std::setprecision(16) <<fParameters[i]<<std::endl;
  }

  void Parametrization::ImportPars(string const& path)
  {
    std::ifstream in(path);
    std::vector<double> inpars;
    double inpar = 0;
    while (in >> inpar)
      inpars.push_back(inpar);

    if (inpars.size() != fNParameters)
    {
      std::cerr << "Error: Not enough parameters in imported vector!" <<std::endl;
      exit(-1);
    };

    for (int i=0; i<fNParameters; i++)
      fParameters[i] = inpars[i];
  }

// ******************** MLP *********************************
/**
 * @brief MultiLayerPerceptron::MultiLayerPerceptron
 * @param settings
 * @param rg
 */
 MultiLayerPerceptron::MultiLayerPerceptron(vector<int> const& arch):
 Parametrization(string("MultiLayerPerceptron")),
 fNLayers(arch.size()),
 fArch(0),
 fWeightMatrix(0),
 fOutputMatrix(0)
 {  
  int  pSize = 0; // size of parameter array
  int* pMap = new int[fNLayers-1]; // map for weight matrix
  
  // Architecture information
  fArch = new int[fNLayers];
  for (int i=0; i<fNLayers; i++)
  {
    fArch[i] = arch[i];
    
    if (i > 0)
    {
      pMap[i-1] = pSize; // startpoint of this layer
      pSize+=fArch[i]*(1+fArch[i-1]); // size of this layer
    }
  }

  // Alloc parameter array
  fParameters = new double[pSize];
  fNParameters = pSize;
  
  // Alloc WeightMatrix (map for parameters)
  fWeightMatrix = new double*[fNLayers - 1];

  // Alloc activation function (output) matrix
  fOutputMatrix = new double*[fNLayers];
  fOutputMatrix[0] = new double[fArch[0]+1];

  for (int i=1; i<fNLayers; i++)
  {
    // point to correct part of parameter array
    fWeightMatrix[i-1] = &fParameters[pMap[i-1]];
    // Alloc each activation function in the layer
    fOutputMatrix[i]   = new double[fArch[i]+1];
  }

 //Init
  fOutputMatrix[0][fArch[0]] =-1.0f;
  for (int i=1; i<fNLayers; i++)
  {
    const int wsz = (fArch[i-1]+1)*fArch[i];
    
    for (int j=0; j<wsz; j++)
      fWeightMatrix[i-1][j] = rng_gaussian(1.0);
    for (int j=0; j<fArch[i]; j++)
      fOutputMatrix[i][j] = 1;
    
    // Threshold term
    fOutputMatrix[i][fArch[i]] = -1.0f;
  }

  
  delete[] pMap;
}

/**
 * @brief MultiLayerPerceptron::MultiLayerPerceptron
 * @param o
 */
 MultiLayerPerceptron::MultiLayerPerceptron(MultiLayerPerceptron const& o):
 Parametrization(o),
 fNLayers(o.fNLayers),
 fArch(0),
 fWeightMatrix(0),
 fOutputMatrix(0)
 {
  int  pSize = 0; // size of parameter array
  int* pMap = new int[fNLayers-1]; // map for weight matrix
  
  // Architecture information
  fArch = new int[fNLayers];
  for (int i=0; i<fNLayers; i++)
  {
    fArch[i] = o.fArch[i];
    
    if (i > 0)
    {
      pMap[i-1] = pSize; // startpoint of this layer
      pSize+=fArch[i]*(1+fArch[i-1]); // size of this layer
    }
  }
  // Alloc WeightMatrix (map for parameters)
  fWeightMatrix = new double*[fNLayers - 1];
  
  // Alloc activation function (output) matrix
  fOutputMatrix = new double*[fNLayers];
  fOutputMatrix[0] = new double[fArch[0]+1];
  
  for (int i=1; i<fNLayers; i++)
  {
    // point to correct part of parameter array
    fWeightMatrix[i-1] = &fParameters[pMap[i-1]];
    // Alloc each activation function in the layer
    fOutputMatrix[i]   = new double[fArch[i]+1];
  }
  
  for (int i=0; i<fNLayers; i++) // Threshold term init
    fOutputMatrix[i][fArch[i]] = -1.0f;
  
  delete[] pMap;
}

/**
 * @brief MultiLayerPerceptron::~MultiLayerPerceptron
 */
 MultiLayerPerceptron::~MultiLayerPerceptron()
 {
  delete[] fOutputMatrix[0];
  for (int i=1; i<fNLayers; i++)
    delete[] fOutputMatrix[i];
  
  delete[] fOutputMatrix;
  
  delete[] fArch;
  delete[] fWeightMatrix;
  delete[] fParameters;
}

/**
 * @brief MultiLayerPerceptron::Duplicate
 * @return
 */
 Parametrization* MultiLayerPerceptron::Duplicate()
 {
  return new MultiLayerPerceptron(*this);
}

/**
 * @brief MultiLayerPerceptron::Compute
 * @param in
 * @param out
 */
 void MultiLayerPerceptron::Compute(const double* in,double* out) const
 {
  // setup input
  for (int i=0; i<fArch[0]; i++)
    fOutputMatrix[0][i] = in[i];
  
  for (int i=1; i<(fNLayers); i++)
    for (int j=0; j<fArch[i]; j++)
    {
      double h=0.0f;
      
      double *p = &fWeightMatrix[i-1][j*(1+fArch[i-1])]; // seems to help the compiler out
      for (int k=0; k<=fArch[i-1]; k++) // <= due to threshold term
        h-= (*(p+k))*fOutputMatrix[i-1][k];
      
      fOutputMatrix[i][j]=1.0f/(1.0f+exp(h));
    }

    for (int i=0; i<fArch[fNLayers-1]; i++)
      out[i] = fOutputMatrix[fNLayers-1][i];
  }

  const int MultiLayerPerceptron::GetNumNodeParams(int const& layer)
  {
    if (layer <=0 || layer >= fNLayers)
    {
      cerr<< "MultiLayerPerceptron::GetNumNodeParams Error: layer requested ("<<layer<<") is out of bounds!"<<endl;
      exit(-1);
    }

    return fArch[layer-1] + 1;
  }

  double* MultiLayerPerceptron::GetNodeParams(int const& layer, int const& node)
  {
    if (layer <=0 || layer >= fNLayers)
    {
      cerr<< "MultiLayerPerceptron::GetNodeParams Error: layer requested ("<<layer<<") is out of bounds!"<<endl;
      exit(-1);
    }

    if (node < 0 || node >= fArch[layer] )
    {
      cerr<< "MultiLayerPerceptron::GetNodeParams Error: node requested ("<<node<<") is out of bounds!"<<endl;
      exit(-1);
    }

    return &fWeightMatrix[layer-1][node];
  }

  void MultiLayerPerceptron::Export(string const& path, string const& inputs) const
  {
    ofstream outputParam(path.c_str());
    for (int i=0; i<fNLayers; i++)
      outputParam << fArch[i]<<"\t";
    
    outputParam<<endl;
    outputParam << inputs<<endl;
    outputParam << "** Node, Previous Node, Threshold, Weight ** "<<endl;

    int iNode = fArch[0];
    int lNode = 0;
    for (int i=1; i<(fNLayers); i++)
    {
      for (int j=0; j<fArch[i]; j++)
      {
        double *p = &fWeightMatrix[i-1][j*(1+fArch[i-1])]; // seems to help the compiler out
        for (int k=0; k<fArch[i-1]; k++) // <= due to threshold term
            outputParam << iNode<<"\t"<<k+lNode<<"\t"<<-*(p+fArch[i-1])<< "\t"<< *(p+k) <<endl; 
        iNode++;
      }

      lNode+=fArch[i-1];
    }

  outputParam.close();
}