#include "TFile.h"
#include "TH1D.h"

void YieldFromHist(TString cat = "pt", TString filename = "output/inclusive/Zg_Unfolding_Subtract_Bkg_data_mu_endcap_isPVGoodstat99_new.root") {

  TFile *infile = new TFile(filename, "read");

  TString hstatname;
  //if (cat.Contains("pt")) hstatname = "histUnfoldStat";
  if (cat.Contains("pt")) hstatname = "histUnfoldStat_ptllg";
  else if (cat.Contains("mllg")) hstatname = "histUnfoldStat_mllg";
  else if (cat.Contains("njet")) hstatname = "histUnfoldStat_njet";

  TString htotname;
  //if (cat.Contains("pt")) htotname = "histUnfoldTotal";
  if (cat.Contains("pt")) htotname = "histUnfoldTotal_ptllg";
  else if (cat.Contains("mllg")) htotname = "histUnfoldTotal_mllg";
  else if (cat.Contains("njet")) htotname = "histUnfoldTotal_njet";


  TH1D *hstat = (TH1D*) infile->Get(hstatname);
  TH1D *htot = (TH1D*) infile->Get(htotname);

  int nbins = hstat->GetNbinsX();

  std::cout << std::setprecision(1) << std::fixed << endl;
  cout << "unfolded yield: " << endl;

  for (int i = 1; i<= nbins; i++) {
    cout << hstat->GetBinContent(i) << " ";
  } 

  cout << endl;
  for (int i = 1; i<= nbins; i++) {
    cout << hstat->GetBinError(i) << " ";
  } 

  cout << "\n error of unfold " << endl;
  for (int i = 1; i<= nbins; i++) {
    float err = sqrt(pow(htot->GetBinError(i),2) - pow(hstat->GetBinError(i),2));
    cout << err << " " ;
  } 

  cout << "\n done!!!!" << endl;



}

