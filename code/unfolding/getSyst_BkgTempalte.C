#include "TFile.h"
#include "TH1D.h"

void getSyst_BkgTempalte(bool ele = 1, bool barrel =1, TString cat = "pt") {


  //for inclusive 
  TString fname_nom = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";

  TString fname_biasUp = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";
  TString fname_biasDn = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";


  if (ele) {
    fname_nom += "_ele";
    fname_biasUp += "_ele";
    fname_biasDn += "_ele";
  }
  else {
    fname_nom += "_mu";
    fname_biasUp += "_mu";
    fname_biasDn += "_mu";
  }

  if (barrel) {
    fname_nom += "_barrel";
    fname_biasUp += "_barrel";
    fname_biasDn += "_barrel";
  }
  else {
    fname_nom += "_endcap";
    fname_biasUp += "_endcap";
    fname_biasDn += "_endcap";
  }

  fname_nom += "_isPVGood_gencut.root";
  cout << "input file: " << fname_nom << endl;
  //template bkg
  //fname_biasUp += "_biasUp_isPVGood_gencut.root";
  //fname_biasDn += "_biasDn_isPVGood_gencut.root";

  fname_biasUp += "_ttg_vv_bkg_UncUp_isPVGood_gencut.root";
  fname_biasDn += "_ttg_vv_bkg_UncDn_isPVGood_gencut.root";

  TFile *fnom = new TFile(fname_nom, "read");
  TFile *fbiasUp = new TFile(fname_biasUp, "read");
  TFile *fbiasDn = new TFile(fname_biasDn, "read");

  TString hname;
  //if (cat.Contains("pt")) hname = "histUnfoldTotal";
  if (cat.Contains("pt")) hname = "histUnfoldTotal_ptllg";
  else if (cat.Contains("mllg")) hname = "histUnfoldTotal_mllg";
  else if (cat.Contains("njet")) hname = "histUnfoldTotal_njet";

  TH1D *hnom = (TH1D*) fnom->Get(hname);
  TH1D *hbiasUp = (TH1D*) fbiasUp->Get(hname);
  TH1D *hbiasDn = (TH1D*) fbiasDn->Get(hname);

  int nbins = hnom->GetNbinsX();

  std::cout << std::setprecision(1) << std::fixed << endl;

  cout << "uncertainty due to background template" << endl;
  for (int i = 1; i<= nbins; i++) {
    float ev_nom = hnom->GetBinContent(i);
    float ev_biasUp = hbiasUp->GetBinContent(i);
    float ev_biasDn = hbiasDn->GetBinContent(i);

    float sys_err = TMath::Max(fabs(ev_nom-ev_biasUp), fabs(ev_nom-ev_biasDn));
    cout << sys_err << " ";
  }
  cout << "\n " << endl;

  gStyle->SetOptStat(0);

}
