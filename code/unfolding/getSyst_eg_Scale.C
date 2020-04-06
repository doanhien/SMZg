#include "TFile.h"
#include "TH1D.h"

void getSyst_eg_Scale(bool ele = 1, bool barrel =1, TString cat = "pt") {

  TString fname_nom = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";
  TString fname_scale = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";
  TString fname_resol = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";

  if (ele) {
    fname_nom += "_ele";
    fname_scale += "_ele";
    fname_resol += "_ele";
  }
  else {
    fname_nom += "_mu";
    fname_scale += "_mu";
    fname_resol += "_mu";
  }

  if (barrel) {
    fname_nom += "_barrel";
    fname_scale += "_barrel";
    fname_resol += "_barrel";
  }
  else {
    fname_nom += "_endcap";
    fname_scale += "_endcap";
    fname_resol += "_endcap";
  }


  fname_nom += "_isPVGood_gencut.root";
  cout << "nominal file " << fname_nom.Data() << endl;

  TFile *fnom = new TFile(fname_nom, "read");

  TFile *fscale_stat_up = new TFile(fname_scale + "_isPVGood_ele_scale_stat_up.root", "read");
  TFile *fscale_stat_dn = new TFile(fname_scale + "_isPVGood_ele_scale_stat_dn.root", "read");
  TFile *fscale_syst_up = new TFile(fname_scale + "_isPVGood_ele_scale_syst_up.root", "read");
  TFile *fscale_syst_dn = new TFile(fname_scale + "_isPVGood_ele_scale_syst_dn.root", "read");
  TFile *fscale_gain_up = new TFile(fname_scale + "_isPVGood_ele_scale_gain_up.root", "read");
  TFile *fscale_gain_dn = new TFile(fname_scale + "_isPVGood_ele_scale_gain_dn.root", "read");

  TString hname;

  //if (cat.Contains("pt")) hname = "histUnfoldOutput";
  if (cat.Contains("pt")) hname = "histUnfoldOutput_ptllg";
  else if (cat.Contains("mllg")) hname = "histUnfoldOutput_mllg";
  else if (cat.Contains("njet")) hname = "histUnfoldOutput_njet";

  TH1D *hnom = (TH1D*) fnom->Get(hname);
  TH1D *hscale_stat_up = (TH1D*) fscale_stat_up->Get(hname);
  TH1D *hscale_stat_dn = (TH1D*) fscale_stat_dn->Get(hname);
  TH1D *hscale_syst_up = (TH1D*) fscale_syst_up->Get(hname);
  TH1D *hscale_syst_dn = (TH1D*) fscale_syst_dn->Get(hname);
  TH1D *hscale_gain_up = (TH1D*) fscale_gain_up->Get(hname);
  TH1D *hscale_gain_dn = (TH1D*) fscale_gain_dn->Get(hname);

  int nbins = hnom->GetNbinsX();

  vector<float> unc_scale;
  vector<float> unc_resol;

  unc_scale.clear();
  unc_resol.clear();

  for (int i = 1; i<= nbins; i++) {
    float ev_nom = hnom->GetBinContent(i);

    float ev_scale[6];
    ev_scale[0] = hscale_stat_up->GetBinContent(i);
    ev_scale[1] = hscale_stat_dn->GetBinContent(i);
    ev_scale[2] = hscale_syst_up->GetBinContent(i);
    ev_scale[3] = hscale_syst_dn->GetBinContent(i);
    ev_scale[4] = hscale_gain_up->GetBinContent(i);
    ev_scale[5] = hscale_gain_dn->GetBinContent(i);

    float max_scale = 0.;
    float max_resol = 0.;

    for (int is = 0; is < 6; is++) {
      if (max_scale < fabs(ev_nom - ev_scale[is])) max_scale = fabs(ev_nom - ev_scale[is]);
    }

    unc_scale.push_back(max_scale);
    unc_resol.push_back(max_resol);

  }

  cout << setprecision(2) << std::fixed << endl;

  for (int i = 1; i <= nbins; i++) {
    cout << hnom->GetBinContent(i) << " ";
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_stat_up->GetBinContent(i) << " " ;
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_stat_dn->GetBinContent(i) << " ";
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_syst_up->GetBinContent(i) << " ";
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_syst_dn->GetBinContent(i) << " ";
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_gain_up->GetBinContent(i) << " ";
  }

  cout << endl;
  for (int i = 1; i <= nbins; i++) {
    cout << hscale_gain_dn->GetBinContent(i) << " ";
  }

  cout << "\n------------uncertainty from scale--------------: " << endl;
  for (int i = 0; i < unc_scale.size(); i++) {
    cout << unc_scale[i] << " ";
  }

  cout << endl;

}
