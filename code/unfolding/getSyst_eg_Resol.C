#include "TFile.h"
#include "TH1D.h"

void getSyst_eg_Resol(bool ele = 1, bool barrel =1, TString cat = "pt") {

  TString fname_nom = "outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_data";
  TString fname_resol = "output_EGM_Resol/Zg_Unfolding_Subtract_Bkg_data";

  if (ele) {
    fname_nom += "_ele";
    fname_resol += "_ele";
  }
  else {
    fname_nom += "_mu";
    fname_resol += "_mu";
  }
  if (barrel) {
    fname_nom += "_barrel";
    fname_resol += "_barrel";
  }
  else {
    fname_nom += "_endcap";
    fname_resol += "_endcap";
  }


  fname_nom += "_isPVGood_gencut.root";

  cout << "nominal file " << fname_nom.Data() << endl;

  TFile *fnom = new TFile(fname_nom, "read");

  TString hname;

  //if (cat.Contains("pt")) hname = "histUnfoldOutput";
  if (cat.Contains("pt")) hname = "histUnfoldOutput_ptllg";
  else if (cat.Contains("mllg")) hname = "histUnfoldOutput_mllg";
  else if (cat.Contains("njet")) hname = "histUnfoldOutput_njet";

  TH1D *hnom = (TH1D*) fnom->Get(hname);
  int nbins = hnom->GetNbinsX();

  ofstream fout;
  string outfname;
  ofstream fout1, fout2, fout3, fout4;
  string outfname1, outfname2, outfname3, outfname4;

  outfname = "EGM_ER_ES_Unc/ele_resol";
  if (cat == "njet") outfname += "_njet";
  else {
    //if (cat == "pt") outfname += "_pt";
    if (cat == "pt") outfname += "_ptllg";
    if (cat == "mllg") outfname += "_mllg";
    if (barrel) outfname += "_barrel";
    else outfname += "_endcap";
  }

  outfname += "_Pho20";

  if (ele) {
    outfname1 = outfname + "_ele_rhoup.txt";
    outfname2 = outfname + "_ele_rhodn.txt";
    outfname3 = outfname + "_ele_phiup.txt";
    outfname4 = outfname + "_ele_phidn.txt";
  }
  else {
    outfname1 = outfname + "_mu_rhoup.txt";
    outfname2 = outfname + "_mu_rhodn.txt";
    outfname3 = outfname + "_mu_phiup.txt";
    outfname4 = outfname + "_mu_phidn.txt";
  }


  if (ele) outfname += "_ele.txt";
  else  outfname += "_mu.txt";

  fout.open(outfname, ios::app);
  fout1.open(outfname1, ios::app);
  fout2.open(outfname2, ios::app);
  fout3.open(outfname3, ios::app);
  fout4.open(outfname4, ios::app);

  Float_t yield[nbins];
  TString ssname = "EGM_ER_ES_Unc/ele_resol_unc";
  //if (cat.Contains("pt")) ssname += "_ptg";
  if (cat.Contains("pt")) ssname += "_ptllg";
  else if (cat.Contains("mllg")) ssname += "_mllg";
  else if  (cat.Contains("njet")) ssname += "_njet";

  if (ele) ssname += "_ele";
  else ssname += "_mu";

  if (barrel) ssname += "_EB";
  else ssname += "_EE";
  ssname += ".root";

  TFile *f = new TFile(ssname,"RECREATE");
  TTree *outtree = new TTree("outtree", "tree for eg resol");
  outtree->Branch("nbins",  &nbins);
  outtree->Branch("yield",  &yield,   "yield[nbins]/F");


  for (int iter = 1; iter <= 100; iter++) {
    TFile *fresol_rho_up = new TFile(fname_resol + Form("_isPVGood_ele_resol_rho_up_%d.root",iter), "read");
    TFile *fresol_rho_dn = new TFile(fname_resol + Form("_isPVGood_ele_resol_rho_dn_%d.root",iter), "read");
    TFile *fresol_phi_up = new TFile(fname_resol + Form("_isPVGood_ele_resol_phi_up_%d.root",iter), "read");
    TFile *fresol_phi_dn = new TFile(fname_resol + Form("_isPVGood_ele_resol_phi_dn_%d.root",iter), "read");


    TH1D *hresol_rho_up = (TH1D*) fresol_rho_up->Get(hname);
    TH1D *hresol_rho_dn = (TH1D*) fresol_rho_dn->Get(hname);
    TH1D *hresol_phi_up = (TH1D*) fresol_phi_up->Get(hname);
    TH1D *hresol_phi_dn = (TH1D*) fresol_phi_dn->Get(hname);

    for (int i = 1; i<= nbins; i++) {
      float ev_nom = hnom->GetBinContent(i);

      float ev_resol[4];
      ev_resol[0] = hresol_rho_up->GetBinContent(i);
      ev_resol[1] = hresol_rho_dn->GetBinContent(i);
      ev_resol[2] = hresol_phi_up->GetBinContent(i);
      ev_resol[3] = hresol_phi_dn->GetBinContent(i);

      fout1 << ev_resol[0] << " ";
      fout2 << ev_resol[1] << " ";
      fout3 << ev_resol[2] << " ";
      fout4 << ev_resol[3] << " ";

      float max_resol = 0.;

      for (int ir = 0; ir <4; ir++) {
        if (max_resol < fabs(ev_nom - ev_resol[ir])) max_resol = fabs(ev_nom - ev_resol[ir]);
      }

      fout << max_resol << " ";
      yield[i-1] = max_resol;

    }

    outtree->Fill();

    fout << endl;
    fout1 << endl;
    fout2 << endl;
    fout3 << endl;
    fout4 << endl;


    fresol_rho_up->Close();
    fresol_rho_dn->Close();
    fresol_phi_up->Close();
    fresol_phi_dn->Close();

  }

  f->Write();

}

