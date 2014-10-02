//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//
// Now the main plotting of histograms routines
// Plot together signal and background


// pt of the HH system
TCanvas* c_pthh = new TCanvas();
vector<TH1D*> histo_pthh;
double const pthh_min=0;
double const pthh_max=500;
int nbin_pthh=20;
double histo_pthh_binwidth = (  pthh_max -  pthh_min ) / nbin_pthh;

// pt of the H candidates
TCanvas* c_pth = new TCanvas();
vector<TH1D*> histo_pth;
double const pth_min=0;
double const pth_max=600;
int nbin_pth=20;
double histo_pth_binwidth = (  pthh_max -  pthh_min ) / nbin_pthh;

// Plot legend
TLegend *leg = new TLegend(0.70,0.68,0.89,0.89,NULL,"brNDC");

void histo_init(){

  // Select here linear or log scales of y axis
  bool ylog=true;
  if(ylog){
    
    // Initialize the various histograms
    
    // pt of the HH system
    c_pthh -> SetLogy();

    // pt of the higgs candidates
    c_pth -> SetLogy();

  }

}

void histo_create(){

  // Create the various histograms
  // This is done for each of the signal or background samples

  std::cout<<"\n ********************************************************************** \n"<<std::endl;

  // pt of the HH system
  histo_pthh.push_back(new TH1D("histo_pthh","histo_pthh",nbin_pthh,  pthh_min,  pthh_max)) ;
   // pt of the h system
  histo_pth.push_back(new TH1D("histo_pth","histo_pth",nbin_pth,  pth_min,  pth_max)) ;

  std::cout<<"\n ********************************************************************** \n"<<std::endl;

}

// The filling is performed in the corresponding part of the code
// Finally the plot the histograms
// Recall that the filling needs to be performed in individual analysis categories

void histo_plot(double scalefactor){

  // Plot of the histograms with corresponding style
  // Save into the various relevant formats
  // Before plotting, we also check that the normalization of the histogram is the correct one
  // Using info on total cross-section and number of events that pass cuts, as done in the
  // original version of the program
  // Also scale by the total number of events
  // scalefactor = 1/nevent
  // In some cases need to renormalize histogram by number of entries??
  
  std::cout<<"\n ********************************************************************** \n"<<std::endl;

  // Here select if we scale to the total cross-section
  // or we want an area of unit 1
  bool histo_normalized=true;
  
  std::cout<<"\n Drawing pthh histogram \n "<<std::endl;
  
  c_pthh->cd();
  histo_pthh.back()->Scale(scalefactor/2); // Two entries per event
  histo_pthh.back()->GetXaxis()->SetTitle("p_{T}^{hh} (GeV) ");
  histo_pthh.back()->GetXaxis()->CenterTitle(true);
  histo_pthh.back()->GetYaxis()->CenterTitle(true);
  histo_pthh.back()->GetYaxis()->SetTitle("d #sigma / d p_{T}^{hh} ( fb / GeV ) ");
  if(histo_normalized) histo_pthh.back()->GetYaxis()->SetTitle("d #sigma / d p_{T}^{hh} ( AU ) ");
  histo_pthh.back()->SetTitle("Gluon Fusion HH, LHC 14 TeV");
  histo_pthh.back()->SetLineWidth(3);
  if(histo_pthh.size()==1)histo_pthh.back()->SetLineColor(1);
  if(histo_pthh.size()==2)histo_pthh.back()->SetLineColor(2);
  if(histo_pthh.size()==3)histo_pthh.back()->SetLineColor(4);
  histo_pthh.back()->SetLineStyle(histo_pthh.size());
  histo_pthh.back()->GetYaxis()->SetRangeUser(1e-4,1.0);
  histo_pthh.back()->GetYaxis()->SetLimits(1e-4,1.0);
  double total_xsec=0;
  for(int j=1; j<=histo_pthh.back()->GetNbinsX(); j++   ){
    total_xsec+= histo_pthh.back()->GetBinContent(j) * histo_pthh_binwidth ;
  }
  std::cout<<"xsec (fb) = "<<total_xsec<<std::endl;
  if(histo_normalized) histo_pthh.back()->Scale(1.0/total_xsec);
  if( histo_pthh.size()==1)histo_pthh.back()->Draw();
  if( histo_pthh.size()>1)histo_pthh.back()->Draw("same");

  std::cout<<"\n Drawing pth histogram \n "<<std::endl;
  
  c_pth->cd();
  histo_pth.back()->Scale(scalefactor/2); // Two entries per event
  histo_pth.back()->GetXaxis()->SetTitle("p_{T}^{h} (GeV) ");
  histo_pth.back()->GetXaxis()->CenterTitle(true);
  histo_pth.back()->GetYaxis()->CenterTitle(true);
  histo_pth.back()->GetYaxis()->SetTitle("d #sigma / d p_{T}^{h} ( fb / GeV ) ");
  if(histo_normalized) histo_pth.back()->GetYaxis()->SetTitle("d #sigma / d p_{T}^{h} ( AU ) ");
  histo_pth.back()->SetTitle("Gluon Fusion HH, LHC 14 TeV");
  histo_pth.back()->SetLineWidth(3);
  if(histo_pth.size()==1)histo_pth.back()->SetLineColor(1);
  if(histo_pth.size()==2)histo_pth.back()->SetLineColor(2);
  if(histo_pth.size()==3)histo_pth.back()->SetLineColor(4);
  histo_pth.back()->SetLineStyle(histo_pth.size());
  histo_pth.back()->GetYaxis()->SetRangeUser(1e-4,1.0);
  histo_pth.back()->GetYaxis()->SetLimits(1e-4,1.0);
  total_xsec=0;
  for(int j=1; j<=histo_pth.back()->GetNbinsX(); j++   ){
    total_xsec+= histo_pth.back()->GetBinContent(j) * histo_pth_binwidth ;
  }
  std::cout<<"xsec (fb) = "<<total_xsec<<std::endl;
  if(histo_normalized) histo_pth.back()->Scale(1.0/total_xsec);
  if( histo_pth.size()==1)histo_pth.back()->Draw();
  if( histo_pth.size()>1)histo_pth.back()->Draw("same");


  std::cout<<"\n ********************************************************************** \n"<<std::endl;

}

void histo_plot_final(){

  // Here we print the histograms into the files with the corresponding format
  // Separately for signal and background events
  
  // Legend
  leg->SetLineStyle(1);
  leg->SetBorderSize(1);
  leg->SetFillColor(0);
  TLegendEntry *entry;

  entry=leg->AddEntry("Graph1","hh SM","L");
  entry->SetLineStyle(1);
  entry->SetLineColor(1);
  entry->SetLineWidth(3);
  entry=leg->AddEntry("Graph1","hh BSM, #lambda = 10*#lambda_{SM} ","L");
  entry->SetLineStyle(2);
  entry->SetLineColor(2);
  entry->SetLineWidth(3);
  entry=leg->AddEntry("Graph1","QCD 4b","L");
  entry->SetLineStyle(3);
  entry->SetLineColor(4);
  entry->SetLineWidth(3);

  // Set the corresponding final state label
  string label="14tev_bbbb";

  c_pthh->cd();
  gStyle->SetOptStat(0);
  leg->Draw();
  string file1="histo_pthh_"+label+".eps";
  string file2="histo_pthh_"+label+".C";
  c_pthh->SaveAs(file1.c_str());
  c_pthh->SaveAs(file2.c_str());

  c_pth->cd();
  gStyle->SetOptStat(0);
  leg->Draw();
  file1="histo_pth_"+label+".eps";
  file2="histo_pth_"+label+".C";
  c_pth->SaveAs(file1.c_str());
  c_pth->SaveAs(file2.c_str());

}

void histo_fill(string histofill, double xsec_fill, double entry_fill){

  // Check that the cross-section is not zero
   if(xsec_fill < 1e-10 ){
    std::cout<<"Stopping, using zero weight in filling histograms"<<std::endl;
    std::cout<<"xsec_fill = "<<xsec_fill<<std::endl;
    std::cout<<"histofill = "<<histofill<<std::endl;
    exit(-10);
   }

   // Check that you are not filling the histogram with zero entries
   // in principle this is not allowed
   if(fabs(entry_fill) < 1e-20 ){
     std::cout<<"Stopping,  filling histograms with zeros, not allowed"<<std::endl;
     std::cout<<"entry_fill = "<<entry_fill<<std::endl;
     std::cout<<"histofill = "<<histofill<<std::endl;
     exit(-10);
   }
   
   if(histofill=="pthh"){
     histo_pthh.back()->Fill(entry_fill,xsec_fill/histo_pthh_binwidth);
   }
   else if(histofill=="pth"){
     histo_pth.back()->Fill(entry_fill,xsec_fill/histo_pth_binwidth);
   }

   else{
     std::cout<<"Invalid option for histogram plotting, histofill = "<<histofill<<std::endl;
     exit(-10);
   }
   
   // End of the routine
   return;

}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
