#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

void unc_Bkg(bool ele = 1, bool eb = 1) {


  TString infname = "TTG_VV_bkg";
  if (eb) infname += "_EBpho";
  else  infname += "_EEpho";
  if (ele) infname += "_eleChan";
  else infname += "_muChan";

  TFile *fbgk = new TFile( infname + ".root", "read");
  TFile *fbgk_up = new TFile(infname + "_xsUp.root", "read");
  TFile *fbgk_dn = new TFile(infname + "_xsDn.root", "read");

  TH1F *hbkg_ttg_vv_pt = (TH1F*) fbgk->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_pt_up = (TH1F*) fbgk_up->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_pt_dn = (TH1F*) fbgk_dn->Get("hbkg_ttg_vv_ptg");

  TH1F *hbkg_ttg_vv_mllg = (TH1F*) fbgk->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_mllg_up = (TH1F*) fbgk_up->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_mllg_dn = (TH1F*) fbgk_dn->Get("hbkg_ttg_vv_mllg");

  TH1F *hbkg_ttg_vv_ptllg = (TH1F*) fbgk->Get("hbkg_ttg_vv_ptllg");
  TH1F *hbkg_ttg_vv_ptllg_up = (TH1F*) fbgk_up->Get("hbkg_ttg_vv_ptllg");
  TH1F *hbkg_ttg_vv_ptllg_dn = (TH1F*) fbgk_dn->Get("hbkg_ttg_vv_ptllg");

  int nbin_pt = hbkg_ttg_vv_pt->GetNbinsX();
  int nbin_mllg = hbkg_ttg_vv_mllg->GetNbinsX();
  int nbin_ptllg = hbkg_ttg_vv_ptllg->GetNbinsX();

  //pho pt
  for (int i = 1; i <= nbin_pt; i++) {
    float diff_up = fabs(hbkg_ttg_vv_pt->GetBinContent(i) - hbkg_ttg_vv_pt_up->GetBinContent(i));
    float diff_dn = fabs(hbkg_ttg_vv_pt->GetBinContent(i) - hbkg_ttg_vv_pt_dn->GetBinContent(i));
    float xs_el = hbkg_ttg_vv_pt->GetBinContent(i) - hbkg_ttg_vv_pt->GetBinError(i) - TMath::Max(diff_up,diff_dn);
    hbkg_ttg_vv_pt->SetBinContent(i,xs_el);
  }

  //mllg
  for (int i = 1; i <= nbin_mllg; i++) {
    float diff_up = fabs(hbkg_ttg_vv_mllg->GetBinContent(i) - hbkg_ttg_vv_mllg_up->GetBinContent(i));
    float diff_dn = fabs(hbkg_ttg_vv_mllg->GetBinContent(i) - hbkg_ttg_vv_mllg_dn->GetBinContent(i));
    float xs_el = hbkg_ttg_vv_mllg->GetBinContent(i) - hbkg_ttg_vv_mllg->GetBinError(i) - TMath::Max(diff_up,diff_dn);
    hbkg_ttg_vv_mllg->SetBinContent(i,xs_el);
  }

  //pt llg
  for (int i = 1; i <= nbin_ptllg; i++) {
    float diff_up = fabs(hbkg_ttg_vv_ptllg->GetBinContent(i) - hbkg_ttg_vv_ptllg_up->GetBinContent(i));
    float diff_dn = fabs(hbkg_ttg_vv_ptllg->GetBinContent(i) - hbkg_ttg_vv_ptllg_dn->GetBinContent(i));
    float xs_el = hbkg_ttg_vv_ptllg->GetBinContent(i) - hbkg_ttg_vv_ptllg->GetBinError(i) - TMath::Max(diff_up,diff_dn);
    hbkg_ttg_vv_ptllg->SetBinContent(i,xs_el);
  }

  
  //write output
  TFile *fout = new TFile(infname + "_uncDn.root", "recreate");
  fout->cd();
  hbkg_ttg_vv_pt->Write("hbkg_ttg_vv_ptg");
  hbkg_ttg_vv_mllg->Write("hbkg_ttg_vv_mllg");
  hbkg_ttg_vv_ptllg->Write("hbkg_ttg_vv_ptllg");


  fout->Write();
  fout->Close();
  
  

  TFile *fbgk_jet = new TFile("TTG_VV_bkg_njet.root", "read");
  TFile *fbgk_jet_up = new TFile("TTG_VV_bkg_njet_xsUp.root", "read");
  TFile *fbgk_jet_dn = new TFile("TTG_VV_bkg_njet_xsDn.root", "read");

  TH1F *hbkg_ttg_vv_njet_ele = (TH1F*) fbgk_jet->Get("hbkg_ttg_vv_njet_ele");
  TH1F *hbkg_ttg_vv_njet_ele_up = (TH1F*) fbgk_jet_up->Get("hbkg_ttg_vv_njet_ele");
  TH1F *hbkg_ttg_vv_njet_ele_dn = (TH1F*) fbgk_jet_dn->Get("hbkg_ttg_vv_njet_ele");

  TH1F *hbkg_ttg_vv_njet_mu = (TH1F*) fbgk_jet->Get("hbkg_ttg_vv_njet_mu");
  TH1F *hbkg_ttg_vv_njet_mu_up = (TH1F*) fbgk_jet_up->Get("hbkg_ttg_vv_njet_mu");
  TH1F *hbkg_ttg_vv_njet_mu_dn = (TH1F*) fbgk_jet_dn->Get("hbkg_ttg_vv_njet_mu");

  int nbin_jet = hbkg_ttg_vv_njet_ele->GetNbinsX();
  for (int i = 1; i <= nbin_jet; i++) {
    float diff_ele_up = fabs(hbkg_ttg_vv_njet_ele->GetBinContent(i) - hbkg_ttg_vv_njet_ele_up->GetBinContent(i));
    float diff_ele_dn = fabs(hbkg_ttg_vv_njet_ele->GetBinContent(i) - hbkg_ttg_vv_njet_ele_dn->GetBinContent(i));
    float xs_el = hbkg_ttg_vv_njet_ele->GetBinContent(i) - hbkg_ttg_vv_njet_ele->GetBinError(i) - TMath::Max(diff_ele_up,diff_ele_dn);
    hbkg_ttg_vv_njet_ele->SetBinContent(i,xs_el);

    float diff_mu_up = fabs(hbkg_ttg_vv_njet_mu->GetBinContent(i) - hbkg_ttg_vv_njet_mu_up->GetBinContent(i));
    float diff_mu_dn = fabs(hbkg_ttg_vv_njet_mu->GetBinContent(i) - hbkg_ttg_vv_njet_mu_dn->GetBinContent(i));
    float xs_mu = hbkg_ttg_vv_njet_mu->GetBinContent(i) - hbkg_ttg_vv_njet_mu->GetBinError(i) - TMath::Max(diff_mu_up,diff_mu_dn);
    hbkg_ttg_vv_njet_mu->SetBinContent(i,xs_mu);

  }


  TFile *outfjet = new TFile("TTG_VV_bkg_njet_uncDn.root", "recreate");
  outfjet->cd();
  hbkg_ttg_vv_njet_ele->Write("hbkg_ttg_vv_njet_ele");
  hbkg_ttg_vv_njet_mu->Write("hbkg_ttg_vv_njet_mu");

  outfjet->Write();
  outfjet->Close();


}
