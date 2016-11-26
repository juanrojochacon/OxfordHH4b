#include "Riostream.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#ifndef ROOT_TList
#include <TList.h>
#endif
#include "TPaveStats.h"

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>

TH1D* ConvertYodaToTH1(string dir, string sample, string file){
  
  ifstream yodafile( (dir+"/"+sample+"/"+file+".yoda").data() );

  if( !yodafile ){
    std::cout << "File "+dir+"/"+sample+"/"+file+".yoda does not exist!" << std::endl;
 return;
  }
  
  std::string line;
  double xlow, xhigh, sumw, sumw2, sumx, sumx2, numEntries;
  vector<double> v_xlow;
  vector<double> v_xhigh;
  vector<double> v_sumw;
  vector<double> v_sumw2;
  vector<double> v_sumx;
  vector<double> v_sumx2;
  vector<double> v_numEntries;
  
  int count = 0;
  bool readData = false;
  while ( !yodafile.eof()) 
  {
    count++;
    
    getline(yodafile,line);
    // if ( count <= 11 ){
    //   //cout << "IGNORE LINE: " << count << " | " << line << std::endl;
    //   continue;
    // }
    //
    // Find line before
    //
    if(line.find("xlow") != std::string::npos){
      readData = true;
      continue;
    }
    if(line.find("END") != std::string::npos){
      //cout << "END" << std::endl;
      break;
    }
    if(readData)
    {
      // cout << line << "\n";
      std::istringstream iss(line);
      iss >> xlow >> xhigh >> sumw >> sumw2 >> sumx >> sumx2 >> numEntries;

      v_xlow.push_back(xlow);
      v_xhigh.push_back(xhigh);
      v_sumw.push_back(sumw);
      v_sumw2.push_back(sumw2);
      v_sumx.push_back(sumx);
      v_sumx2.push_back(sumx2);
      v_numEntries.push_back(numEntries);
    }
  }
  
  // for (int i=0; i < v_xlow.size(); i++)
  //   cout << v_xlow.at(i) << " | " << v_xhigh.at(i) << " | "<< v_xlow.size() << endl;
  
  int nbins   = 0; 
  nbins   = v_xlow.size();
  int loedge  = 0; 
  loedge  = v_xlow[0];
  int hiedge  = 0; 
  hiedge  = v_xhigh[v_xhigh.size() - 1];
  
  // std::cout << loedge << std::endl;
  
  TH1D* h1 = new TH1D((sample+"_"+file).data(),"h",nbins,loedge,hiedge);
  for (int i=0; i < v_xlow.size(); i++){
    h1->SetBinContent(i+1,v_sumw.at(i));
    h1->SetBinError(i+1,v_sumx.at(i));
  }

  return h1;   
}
void SetColorTH1D(TH1D* h, int col){
  
  h->SetLineColor(1);
  h->SetLineWidth(2);
  h->SetMarkerColor(col);
  h->SetFillColor(col);

}
  
void FIKRI_MakeTH1(){
  
  // typedef map< string, map< string, map< string, map< string, TH1D* > > > > histograms;
  // histograms[btag]
  // int nAnalysis = 3; int nCuts = 4; int nTagCats = 2;
  // string aString[nAnalysis] = {"_res", "_inter", "_boost"};
  // string cString[nCuts] = {"_RCO", "_SIG", "_SDBA", "_SDBB"};
  // string bTagString[nTagCats] = {"_2tag", "_4tag"};
  
  gStyle->SetOptStat(0);
  gROOT->Reset();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
  gStyle->SetPadLeftMargin(0.15);

  // string dir = "./baseline_noPU_atlas_qcd";
  // string file = "histo_m_HH_res_SIG_2tag";
  
  // TH1D* h_4b_QCD_4b   = ConvertYodaToTH1(dir,"SHERPA_QCD4b",file);
  // TH1D* h_4b_QCD_2b2j = ConvertYodaToTH1(dir,"SHERPA_QCD2b2j",file);
  // TH1D* h_4b_QCD_4j   = ConvertYodaToTH1(dir,"SHERPA_QCD4j",file);
  // TH1D* h_4b_HH       = ConvertYodaToTH1(dir,"diHiggs",file);

  // SetColorTH1D(h_4b_QCD_4b   ,2);
  // SetColorTH1D(h_4b_QCD_2b2j ,3);
  // SetColorTH1D(h_4b_QCD_4j   ,4);
  // SetColorTH1D(h_4b_HH       ,0);
  
  // double xLat = 0.2;
  // double yLat = 0.92;
  // double xLeg = xLat + 0.35;
  // double yLeg = yLat;
              
  // double leg_h =  0.03 * 4;
  // TLegend* leg = new TLegend( xLeg, yLeg - leg_h, xLeg + 0.35, yLeg );
  // //leg->SetNColumns( 2 );
  // leg->SetFillStyle( 0 );
  // leg->SetBorderSize( 0 );
  // leg->SetTextFont( 43 );
  // leg->SetTextSize( 18 );
  
  // THStack *bkgd = new THStack("bkgd","");  
  // THStack *sig = new THStack("sig","");
  
  // bkgd->Add(h_4b_QCD_4b);  leg->AddEntry(h_4b_QCD_4b,  "QCD_4b"  , "f");
  // bkgd->Add(h_4b_QCD_2b2j);leg->AddEntry(h_4b_QCD_2b2j,"QCD_2b2j", "f");
  // bkgd->Add(h_4b_QCD_4j);  leg->AddEntry(h_4b_QCD_4j,  "QCD_4j"  , "f");
  
  // h_4b_HH->Scale(5000);     
  // sig->Add(h_4b_HH);       
  // leg->AddEntry(h_4b_HH,  "HH"  , "l");
  
  // TCanvas* c1 = new TCanvas("c","c",600,600);
  // bkgd->Draw("hist");
  // sig->Draw("nostack hist same");
  // leg->Draw();
  
  // TString pdfname = file.data();
  // c1->Print(pdfname+".pdf");
  
  return;
}

