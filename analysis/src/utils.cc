#include "utils.h"
#include "settings.h"

#include <vector>

using namespace std;

 double getDPhi(double phi1, double phi2){
   
   const double PI = 3.14159265359;
   double deltaPhi = phi1-phi2;
   if (deltaPhi > PI)  deltaPhi = (deltaPhi-2*PI);
   if (deltaPhi < -PI) deltaPhi = (deltaPhi+2*PI);
   return deltaPhi;
 }

// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain splitting scales
double SplittingScales( fastjet::PseudoJet const& jet )
{
  if (!jet.has_constituents())
  {
    std::cerr << "ERROR! Splittings (d12, d23, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
    exit(-1);
  }
   
  vector<fastjet::PseudoJet> const& constits = jet.constituents();
      
  fastjet::JetDefinition ekt_jd = fastjet::JetDefinition( fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best);
  const fastjet::ClusterSequence kt_seq_excl = fastjet::ClusterSequence( constits, ekt_jd);
  fastjet::PseudoJet kt_jet = sorted_by_pt( kt_seq_excl.inclusive_jets())[0];
  
  const double split12 = 1.5*sqrt( kt_seq_excl.exclusive_subdmerge( kt_jet, 1));
  return split12; 
}

std::vector< double > SplittingScales( std::vector<fastjet::PseudoJet> const& jetVec )
{
  //vectors that contain the respective splitting scales for all jets
   std::vector<double> split12_vec;
   
   for( int i = 0; i < (int) jetVec.size(); i++)
      split12_vec.push_back(SplittingScales(jetVec[i]));
   
   return split12_vec;
}

// ----------------------------------------------------------------------------------
//Code snippet modified from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetSubstructureVariables#N_subjettiness
double NSubjettiness( fastjet::PseudoJet const& jet, double const& jet_rad )
{  
  double alpha=1;
  if (!jet.has_constituents())
  {
    std::cerr << "ERROR! Subjettiness can only be calculated on jets for which the constituents are known."<< std::endl;
    exit(-1);
  }

  vector<fastjet::PseudoJet> const& constits = jet.constituents();
  if(constits.size()==0)
  {
     std::cerr << "ERROR! Empty jet in Nsubjettiness!"<< std::endl;
     exit(-1);
  }
      
   const double R_kt = 1.5; // ******************************************  SHOULD THIS BE jet_rad? ********************************
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,R_kt,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence kt_clust_seq(constits, jet_def);

   double tau1 = 0; // 1-subjettiness 
   if (constits.size() > 0) // scoping
   {
     double tauNum = 0.0;
     double tauDen = 0.0;

     vector<fastjet::PseudoJet> kt1axes = kt_clust_seq.exclusive_jets(1);
     for (int i = 0; i < (int)constits.size(); i++) 
     {
        // find minimum distance
        double minR = std::numeric_limits<double>::infinity();
        for (int j = 0; j < (int)kt1axes.size(); j++) 
        {
           double tempR = sqrt(constits[i].squared_distance(kt1axes[j])); // delta R distance
           if (tempR < minR) minR = tempR;
        }
        tauNum += constits[i].perp() * pow(minR,alpha);
        tauDen += constits[i].perp() * pow(jet_rad,alpha);
     }
     tau1 = tauNum/tauDen;
    }

   double tau2 = 0; // 2-subjettiness
   if(constits.size()>1)
   {
      double tauNum = 0.0;
      double tauDen = 0.0;

      vector<fastjet::PseudoJet> kt2axes = kt_clust_seq.exclusive_jets(2);
      for (int i = 0; i < (int)constits.size(); i++) 
      {
        // find minimum distance
        double minR = std::numeric_limits<double>::infinity(); // large number
        for (int j = 0; j < (int)kt2axes.size(); j++) 
        {
           double tempR = sqrt(constits[i].squared_distance(kt2axes[j]));
           if (tempR < minR) minR = tempR;
        }
        tauNum += constits[i].perp() * pow(minR,alpha);
        tauDen += constits[i].perp() * pow(jet_rad,alpha);
      }
      tau2 = tauNum/tauDen;
   }
      
  return tau2/tau1;
}

std::vector< double > NSubjettiness( std::vector<fastjet::PseudoJet> const& jetVec, double const& jet_rad )
{
   std::vector<double> tau21_vec;
   for( size_t i = 0; i < std::min(jetVec.size(), (size_t)2); i++)
    tau21_vec.push_back( NSubjettiness(jetVec[i], jet_rad) );
     
   return tau21_vec;
}


// ************************************ Energy correlations *******************************************

double ECF(size_t const& N, double const& beta, fastjet::PseudoJet const& jet)
{
  if (!jet.has_constituents())
  {
     std::cerr << "ERROR! ECF can only be calculated on jets for which the constituents are known."<< std::endl;
     exit(-1);
  }
  
  vector<fastjet::PseudoJet> const& constituents = jet.constituents();
  if(constituents.size()<N){
    return 0;
  }

  double ECF = 0.0;

  switch (N)
  {
    case 0:
      ECF = 1;
      break;

    case 1:
      for (size_t i=0; i<constituents.size(); i++)
        ECF += constituents[i].pt();
      break;

    case 2:
      for (size_t i=0; i<constituents.size(); i++)
        for (size_t j=0; j<i; j++)
        {
          const double dR = std::pow(constituents[i].delta_R(constituents[j]),beta);
          ECF += constituents[i].pt()*constituents[j].pt()*dR;
        }
      break;

    case 3:
      for (size_t i=0; i<constituents.size(); i++)
        for (size_t j=0; j<i; j++)
          for (size_t k=0; k<j; k++)
          {
            const double dR_ij = constituents[i].delta_R(constituents[j]);
            const double dR_ik = constituents[i].delta_R(constituents[k]);
            const double dR_jk = constituents[j].delta_R(constituents[k]);
            const double dR = std::pow(dR_ij*dR_ik*dR_jk, beta);
            ECF += constituents[i].pt()*constituents[j].pt()*constituents[k].pt()*dR;
          }
      break;

    default:
     std::cerr << "ECF::ERROR! Correlation function for N="<<N<<" not implemented!" << std::endl;
     exit(-1);

  }

  return ECF;
}

double LST_C2(double const& beta, fastjet::PseudoJet const& jet)
{
  const double ECF1 = ECF(1, beta, jet);
  const double ECF2 = ECF(2, beta, jet);
  const double ECF3 = ECF(3, beta, jet);

  return (ECF1*ECF3)/(ECF2*ECF2);
}


// ----------------------------------------------------------------------------------
// Recluster with kt algorithm to obtain nsubjettiness
std::vector< double > NSubjettiness( std::vector<fastjet::PseudoJet> const& jetVec, double const& jet_Rmax, double const& jet_Rmin, double const& jet_Rho )
{

   // Reclustering with VR kt algorithm to obtain nsubjettiness
   
   //vector that contain the respective nsubjettiness variables for all jets
   std::vector<double> tau1_vec;
   std::vector<double> tau2_vec;
   std::vector<double> tau3_vec;
   
   std::vector<double> tau21_vec;
   
   double alpha=1;
   
   for( int i = 0; i < (int) jetVec.size(); i++){
   
      // For now: Calculate substructure information only for the two leading jets
      if( i > 1 ) continue;
   
      double tau1 = -1.;
      double tau2 = -1.;
      double tau3 = -1.;
      
      // Calculate effective jet radius
      double jet_Pt = jetVec[i].pt();
      
      double jet_rad;
      if( jet_Pt > jet_Rmax ) jet_rad = jet_Rmax;
      else if( jet_Pt < jet_Rmin ) jet_rad = jet_Rmin;
      else jet_rad = jet_Rho / jet_Pt;

      if (!jetVec.at(i).has_constituents()){
         std::cout << "ERROR! NSubjettiness (tau1, tau2, ...) can only be calculated on jets for which the constituents are known."<< std::endl;
         
         tau1_vec.push_back(-1);
   tau2_vec.push_back(-1);
   tau3_vec.push_back(-1);
  
   //double tau21 = tau2/tau1;
   tau21_vec.push_back(-1);
         continue;
      }
      
      //Code snippet taken from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetSubstructureVariables#N_subjettiness
      vector<fastjet::PseudoJet> constits = jetVec.at(i).constituents();
      if(constits.size()>0){
      
   double R_kt = 1.5; // ******************************************  SHOULD THIS BE jet_rad? ********************************
   fastjet::JetDefinition jet_def = fastjet::JetDefinition(fastjet::kt_algorithm,R_kt,fastjet::E_scheme,fastjet::Best);
   fastjet::ClusterSequence kt_clust_seq(constits, jet_def);
   vector<fastjet::PseudoJet> kt1axes = kt_clust_seq.exclusive_jets(1);
   double tauNum = 0.0;
   double tauDen = 0.0;
   for (int i = 0; i < (int)constits.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (int j = 0; j < (int)kt1axes.size(); j++) {
         double tempR = sqrt(constits[i].squared_distance(kt1axes[j])); // delta R distance
         if (tempR < minR) minR = tempR;
      }
      tauNum += constits[i].perp() * pow(minR,alpha);
      tauDen += constits[i].perp() * pow(jet_rad,alpha);
   }
   tau1 = tauNum/tauDen;

   if(constits.size()>1){
      vector<fastjet::PseudoJet> kt2axes = kt_clust_seq.exclusive_jets(2);
      tauNum = 0.0;
      tauDen = 0.0;
      for (int i = 0; i < (int)constits.size(); i++) {
      // find minimum distance
      double minR = 10000.0; // large number
      for (int j = 0; j < (int)kt2axes.size(); j++) {
         double tempR = sqrt(constits[i].squared_distance(kt2axes[j]));
         if (tempR < minR) minR = tempR;
      }
      tauNum += constits[i].perp() * pow(minR,alpha);
      tauDen += constits[i].perp() * pow(jet_rad,alpha);
      }
      tau2 = tauNum/tauDen;
      
      if(constits.size() > 2){
      
         vector<fastjet::PseudoJet> kt3axes = kt_clust_seq.exclusive_jets(3);
         tauNum = 0.0;
         tauDen = 0.0;
         for (int i = 0; i < (int)constits.size(); i++) {
         // find minimum distance
         double minR = 10000.0; // large number
         for (int j = 0; j < (int)kt3axes.size(); j++) {
      double tempR = sqrt(constits[i].squared_distance(kt3axes[j]));
      if (tempR < minR) minR = tempR;
         }
         tauNum += constits[i].perp() * pow(minR,alpha);
         tauDen += constits[i].perp() * pow(jet_rad,alpha);
         }
         tau3 = tauNum/tauDen;
      }

   }
      }
      
      tau1_vec.push_back(tau1);
      tau2_vec.push_back(tau2);
      tau3_vec.push_back(tau3);
      
      double tau21 = tau2/tau1;
      tau21_vec.push_back(tau21);
   }
   
   return tau21_vec;
}

// ---------------------------------------------------------------------------------------
// Associate charged small-R ("track") jets to large-R ("calo") jet
void get_assoc_trkjets( fastjet::PseudoJet calojet, std::vector<fastjet::PseudoJet> trkjets, std::vector<fastjet::PseudoJet> &matched_trkjets, bool debug=false){

    // vector to hold input clusters and ghosts
    std::vector<fastjet::PseudoJet> input_particles;
    input_particles.clear();

    // jet clusters from large-R jet
    vector<fastjet::PseudoJet> constituents = calojet.constituents();

    if( debug ) std::cout << "calo constituents size = " << constituents.size() << std::endl;

    for(unsigned int constItr = 0; constItr < constituents.size(); ++constItr) {
    
      fastjet::PseudoJet noghost = constituents.at(constItr);
      
      noghost.reset_PtYPhiM (noghost.pt(), noghost.rapidity(), noghost.phi(), 0.0);
      if(noghost.E()<0.) continue;

      // set user index for calo clusters to -1 to differentiate them from "track" constituents later
      noghost.set_user_index(-1);
      input_particles.push_back(noghost);
    }

    if( debug ) std::cout << "calo only input particles size = " << input_particles.size() << std::endl;

    // make ghost PseudoJets out of track jet direction
    for(unsigned int trackJetItr = 0; trackJetItr < trkjets.size(); ++trackJetItr){
    
      fastjet::PseudoJet myghost = trkjets.at(trackJetItr);  
      //if( myghost.pt() <= 20.0 || fabs( myghost.rapidity() ) >= 2.5 ) continue;
      
      myghost.reset_PtYPhiM (1e-12, myghost.rapidity(), myghost.phi(), 0.0);
      if(myghost.E()<0.) continue;

      myghost.set_user_index(trackJetItr);
      input_particles.push_back(myghost);
    }

    if( debug ) std::cout << "calo+track jets input particles size = " << input_particles.size() << std::endl;

    // do ghost association and get list of pseudojet track jets that are associated
    double Rparam = 1.0;
    fastjet::Strategy strategy = fastjet::Best;  // according to atlas reco
    fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;  // according to atlas reco
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, Rparam, recomb_scheme, strategy);

    // run the jet clustering with the above jet definition
    fastjet::ClusterSequence clust_seq(input_particles, jet_def);
    vector<fastjet::PseudoJet> sorted_jets = fastjet::sorted_by_pt( clust_seq.inclusive_jets() );

    if( debug ) std::cout << "number of sorted jets = " << sorted_jets.size() << std::endl;

    fastjet::PseudoJet newJet = sorted_jets.at(0); //there are more jets in the vector, but they all have pT ~0
    if( debug ) std::cout << "new jet constituent size = " << newJet.constituents().size() << std::endl;
    vector<fastjet::PseudoJet> newJet_constituents = newJet.constituents();
    for(unsigned int i=0; i<newJet_constituents.size(); ++i){
      fastjet::PseudoJet constit = newJet_constituents.at(i);
      if( debug ) std::cout << " user index = " << constit.user_index() << ", pt of constit = " << constit.pt() << std::endl;
      if(constit.user_index()>=0){
	  int iter = constit.user_index();
// 	  if(trkjets.at(iter).pt() > 20. && fabs(trkjets.at(iter).eta()) < 2.5 ) matched_trkjets.push_back(trkjets.at(iter));
	  matched_trkjets.push_back(trkjets.at(iter));
      }
    }

    return;
}

// ---------------------------------------------------------------------------------------
// B-tagging based on ATLAS PubNotes 
// (ATL-PHYS-PUB-2014-013, ATL-PHYS-PUB-2014-014)

double btag_eff( double jet_pt ){

  if( jet_pt <= 20. ) return 0.05;
  else if( jet_pt > 20. && jet_pt <= 50. ) return 0.35;
  else if( jet_pt > 50. && jet_pt <= 70. ) return 0.65;
  else if( jet_pt > 70. && jet_pt <= 100. ) return 0.70;
  else if( jet_pt > 100. && jet_pt <= 200. ) return 0.65;
  else if( jet_pt > 200. && jet_pt <= 350. ) return 0.60;
  else if( jet_pt > 350. && jet_pt <= 500. ) return 0.55;
  else return 0.50;
}

double mistag_eff( double jet_pt ){

  if( jet_pt <= 70. ) return 1./150.;
  else if( jet_pt > 70. && jet_pt <= 180. ) return 1./170.;
  else if( jet_pt > 180. && jet_pt <= 250. ) return 1./130.;
  else return 1./100.;
}

double charm_eff( double jet_pt ){

  if( jet_pt <= 70. ) return 1./6.;
  else if( jet_pt > 70. && jet_pt <= 180. ) return 1./4.;
  else if( jet_pt > 180. && jet_pt <= 250. ) return 1./4.5;
  else return 1./5.;
}

static float ranf()/* ranf() is uniform in 0..1 */
{
  return  (float) rand() / (double) RAND_MAX; 
};

// Gaussian random numbers
float box_muller(float m, float s)  /* normal random variate generator */
{               /* mean m, standard deviation s */
  float x1, x2, w, y1;
  static float y2;
  static int use_last = 0;

  if (use_last)           /* use value from previous call */
  {
    y1 = y2;
    use_last = 0;
  }
  else
  {
    do {
      x1 = 2.0 * ranf() - 1.0;
      x2 = 2.0 * ranf() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );

    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return( m + y1 * s );
}

