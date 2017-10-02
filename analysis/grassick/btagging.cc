#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace fastjet;
using namespace std;
/********************************************************************************************************

THIS IS THE UPDATED SCRIPT FOR CALCULATING TOTAL B-TAGGING PROBABILITIES AND NEEDS TO BE INSERTED INTO
THE CODE AND THEN ALL SUBSEQUENT USES OF THE FUNCTION "btagProb" MUST BE ALTERED
AS THE PARAMETERS IT TAKES ARE NOW DIFFERENT.

********************************************************************************************************/

//This is the probability of correctly tagging a b jet as a b
double btag_prob(fastjet::PseudoJet jet){
  double pt = jet.pt();
  if (pt <= 470){// 10th-order polynomial fit
    return -2.9820*pow(10,-24)*pow(pt,10) + 8.0061*pow(10,-21)*pow(pt,9) - 9.3539*pow(10,-18)*pow(pt,8)
      + 6.2359*pow(10,-15)*pow(pt,7) - 2.6140*pow(10,-12)*pow(pt,6) + 7.1605*pow(10,-10)*pow(pt,5)
      - 1.2914*pow(10,-7)*pow(pt,4) + 1.5088*pow(10,-5)*pow(pt,3) - 1.0953*pow(10,-3)*pow(pt,2)
      + 4.5102*pow(10,-2)*pt -5.5887*pow(10,-3);}
  else {//exponential fit
    return exp(-0.0012675*pt + 0.0842395661);}
}



// This is the probability of mistagging a light jet as a b
double btag_mistag(fastjet::PseudoJet jet){
  double pt = jet.pt();
  if (pt <= 300){// 10th-order polynomial fit
    return -1.2957*pow(10,-26)*pow(pt,10) + 4.5327*pow(10,-23)*pow(pt,9) -6.7665*pow(10,-20)*pow(pt,8)
      +5.6403*pow(10,-17)*pow(pt,7) - 2.8639*pow(10,-14)*pow(pt,6) + 9.1997*pow(10,-12)*pow(pt,5)
      - 1.8661*pow(10,-9)*pow(pt,4) + 2.3237*pow(10,-7)*pow(pt,3) - 1.6536*pow(10,-5)*pow(pt,2)
      +5.7809*pow(10,-4)*pt + 7.8167*pow(10,-6);}
  else {
    return 0.0109717;}
}

// This is the probability of mistagging a c jet as a b
double ctag_prob(fastjet::PseudoJet jet){
  double pt = jet.pt();
  if (pt <= 435.334){
    return 5.64989*pow(10,-26)*pow(pt,10) - 8.30785*pow(10,-23)*pow(pt,9) - 2.45638*pow(10,-20)*pow(pt,8)
      + 1.31793*pow(10,-16)*pow(pt,7) - 1.20918*pow(10,-13)*pow(pt,6) + 5.66046*pow(10,-11)*pow(pt,5)
      - 1.55302*pow(10,-8)*pow(pt,4) + 2.55755*pow(10,-6)*pow(pt,3) - 2.44367*pow(10,-4)*pow(pt,2)
      + 1.20369*pow(10,-2)*pt + 2.17015*pow(10,-4);}
     
  else {
    return exp(-0.0020757*pt -1.10695);}
 
}


// function returning all subsets of size k of a vector of integers "vec" stored in a 
// vector of vectors called "result"
 vector<int> v;
//idx is the index it starts at, always use 0
// result is the variable passed as a reference, which is changed by the function
void subset(vector<int> vec,int k,int idx, vector<vector<int> > &result){
    int n = vec.size();
   if(idx==n)
 return;

if(k==1){
    for(int i=idx;i<n;i++)
     {
        v.push_back(vec[i]);
        result.push_back(v);
        v.pop_back();
     }
}

 for(int j=idx;j<n;j++) {
  v.push_back(vec[j]);
  subset(vec,k-1,j+1,result);
  v.pop_back();
  }
 }


//added the last two arguments, and removed nB, nC and nL
//so need to alter every usage of this fn in the existing code accordingly

double btagProb(int const& nTag, vector<PseudoJet> jets, vector<btagType> btags) {

  //creating a vector containing the probabilities that each jet is tagged as a b  
  vector<double> bprob_vec;
  //creating a vector containing the jet indices
  vector<int> indices;

   if (jets.size() != btags.size()){
      return 1;
    }
  for (unsigned int i = 0; i < jets.size(); i++){
    if (btags[i] == BTAG){
      bprob_vec.push_back( btag_prob(jets[i]) );
    }
    else if (btags[i] == CTAG){
      bprob_vec.push_back( ctag_prob(jets[i]) );
    }
    else {
      bprob_vec.push_back( btag_mistag(jets[i]) );
    }
    indices.push_back(i);
    
  }
  
  vector<vector<int> > bperm;
  subset(indices, nTag, 0, bperm);


  vector<double> bprob, nbprob;
  for (unsigned int j = 0; j < bperm.size(); j++){
    bprob.push_back(1);
    nbprob.push_back(1);
    
    for (unsigned int k = 0; k < indices.size(); k++){

      // if the index is not in b_perm[j] we multiply into our vector of 
      // not-btag probabilities
      if ( std::find(bperm[j].begin(), bperm[j].end(),k) == bperm[j].end() ){
	nbprob[j] *= ( 1 - bprob_vec[k]);
      }

      // if the index is in b_perm[j] we multiply into our vector of 
      // btag probabilities
      else {
	bprob[j] *= bprob_vec[k];
      }
    }
  }
    
    double total_prob = 0;

    for (unsigned int r = 0; r < bprob.size(); r++){
      total_prob += ( bprob[r] * nbprob[r] );
         }

    return total_prob;
}
  
	

