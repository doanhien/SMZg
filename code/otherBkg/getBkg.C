#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <iostream>

using namespace std;

void getBkg(bool isEle = 1, bool barrel = 1, int njet = 0) {

  TFile *fttg = new TFile("../ana_jet/minitrees/TTGamma_NoPtCut.root", "read");
  TFile *fWWToLNuQQ = new TFile("../ana_jet/minitrees/WWToLNuQQ_NoPtCut.root", "read");
  TFile *fWWTo2L2Nu = new TFile("../ana_jet/minitrees/WWTo2L2Nu_NoPtCut.root", "read");
  TFile *fWZTo3LNu = new TFile("../ana_jet/minitrees/WZTo3LNu_NoPtCut.root", "read");
  TFile *fWZTo2L2Q = new TFile("../ana_jet/minitrees/WZTo2L2Q_NoPtCut.root", "read");
  TFile *fZZTo2L2Nu = new TFile("../ana_jet/minitrees/ZZTo2L2Nu_NoPtCut.root", "read");
  TFile *fZZTo2L2Q = new TFile("../ana_jet/minitrees/ZZTo2L2Q_NoPtCut.root", "read");
  TFile *fZZTo4L = new TFile("../ana_jet/minitrees/ZZTo4L_NoPtCut.root", "read");

  TTree *ttg = (TTree*)fttg->Get("outtree");
  TTree *tWWToLNuQQ = (TTree*) fWWToLNuQQ->Get("outtree");
  TTree *tWWTo2L2Nu = (TTree*) fWWTo2L2Nu->Get("outtree");
  TTree *tWZTo3LNu = (TTree*) fWZTo3LNu->Get("outtree");
  TTree *tWZTo2L2Q = (TTree*) fWZTo2L2Q->Get("outtree");
  TTree *tZZTo2L2Nu = (TTree*) fZZTo2L2Nu->Get("outtree");
  TTree *tZZTo2L2Q = (TTree*) fZZTo2L2Q->Get("outtree");
  TTree *tZZTo4L = (TTree*) fZZTo4L->Get("outtree");

  const int nbin_pt = 22;
  float pt[nbin_pt+1] = {20, 22, 25, 27, 29, 32, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 105, 120, 200, 1000};

  TH1F *hbkg_ttg_pt = new TH1F("hbkg_ttg_pt", "hbkg_ttg_pt", nbin_pt, pt);
  TH1F *hbkg_WWToLNuQQ_pt = new TH1F("hbkg_WWToLNuQQ_pt", "hbkg_WWToLNuQQ_pt", nbin_pt, pt);
  TH1F *hbkg_WWTo2L2Nu_pt = new TH1F("hbkg_WWTo2L2Nu_pt", "hbkg_WWTo2L2Nu_pt", nbin_pt, pt);
  TH1F *hbkg_WZTo3LNu_pt = new TH1F("hbkg_WZTo3LNu_pt", "hbkg_WZTo3LNu_pt", nbin_pt, pt);
  TH1F *hbkg_WZTo2L2Q_pt = new TH1F("hbkg_WZTo2L2Q_pt", "hbkg_WZTo2L2Q_pt", nbin_pt, pt);
  TH1F *hbkg_ZZTo2L2Nu_pt = new TH1F("hbkg_ZZTo2L2Nu_pt", "hbkg_ZZTo2L2Nu_pt", nbin_pt, pt);
  TH1F *hbkg_ZZTo2L2Q_pt = new TH1F("hbkg_ZZTo2L2Q_pt", "hbkg_ZZTo2L2Q_pt", nbin_pt, pt);
  TH1F *hbkg_ZZTo4L_pt = new TH1F("hbkg_ZZTo4L_pt", "hbkg_ZZTo4L_pt", nbin_pt, pt);
  
  const int nbin_mllg = 22;
  float mllg[nbin_mllg+1] = {70, 85, 88, 92, 95, 100, 105, 110, 115, 120, 135, 150, 170, 190, 210, 240, 270, 300, 350, 400, 470, 640, 3000};

  TH1F *hbkg_ttg_mllg = new TH1F("hbkg_ttg_mllg", "hbkg_ttg_mllg", nbin_mllg, mllg);
  TH1F *hbkg_WWToLNuQQ_mllg = new TH1F("hbkg_WWToLNuQQ_mllg", "hbkg_WWToLNuQQ_mllg", nbin_mllg, mllg);
  TH1F *hbkg_WWTo2L2Nu_mllg = new TH1F("hbkg_WWTo2L2Nu_mllg", "hbkg_WWTo2L2Nu_mllg", nbin_mllg, mllg);
  TH1F *hbkg_WZTo3LNu_mllg = new TH1F("hbkg_WZTo3LNu_mllg", "hbkg_WZTo3LNu_mllg", nbin_mllg, mllg);
  TH1F *hbkg_WZTo2L2Q_mllg = new TH1F("hbkg_WZTo2L2Q_mllg", "hbkg_WZTo2L2Q_mllg", nbin_mllg, mllg);
  TH1F *hbkg_ZZTo2L2Nu_mllg = new TH1F("hbkg_ZZTo2L2Nu_mllg", "hbkg_ZZTo2L2Nu_mllg", nbin_mllg, mllg);
  TH1F *hbkg_ZZTo2L2Q_mllg = new TH1F("hbkg_ZZTo2L2Q_mllg", "hbkg_ZZTo2L2Q_mllg", nbin_mllg, mllg);
  TH1F *hbkg_ZZTo4L_mllg = new TH1F("hbkg_ZZTo4L_mllg", "hbkg_ZZTo4L_mllg", nbin_mllg, mllg);

  const int nbin_ptllg = 16;
  float ptllg[nbin_ptllg+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 90, 120, 150, 1000};

  TH1F *hbkg_ttg_ptllg = new TH1F("hbkg_ttg_ptllg", "hbkg_ttg_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_WWToLNuQQ_ptllg = new TH1F("hbkg_WWToLNuQQ_ptllg", "hbkg_WWToLNuQQ_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_WWTo2L2Nu_ptllg = new TH1F("hbkg_WWTo2L2Nu_ptllg", "hbkg_WWTo2L2Nu_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_WZTo3LNu_ptllg = new TH1F("hbkg_WZTo3LNu_ptllg", "hbkg_WZTo3LNu_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_WZTo2L2Q_ptllg = new TH1F("hbkg_WZTo2L2Q_ptllg", "hbkg_WZTo2L2Q_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_ZZTo2L2Nu_ptllg = new TH1F("hbkg_ZZTo2L2Nu_ptllg", "hbkg_ZZTo2L2Nu_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_ZZTo2L2Q_ptllg = new TH1F("hbkg_ZZTo2L2Q_ptllg", "hbkg_ZZTo2L2Q_ptllg", nbin_ptllg, ptllg);
  TH1F *hbkg_ZZTo4L_ptllg = new TH1F("hbkg_ZZTo4L_ptllg", "hbkg_ZZTo4L_ptllg", nbin_ptllg, ptllg);


  hbkg_ttg_pt->Sumw2();
  hbkg_WWToLNuQQ_pt->Sumw2();
  hbkg_WWTo2L2Nu_pt->Sumw2();
  hbkg_WZTo3LNu_pt->Sumw2();
  hbkg_WZTo2L2Q_pt->Sumw2();
  hbkg_ZZTo2L2Nu_pt->Sumw2();
  hbkg_ZZTo2L2Q_pt->Sumw2();
  hbkg_ZZTo4L_pt->Sumw2();

  hbkg_ttg_mllg->Sumw2();
  hbkg_WWToLNuQQ_mllg->Sumw2();
  hbkg_WWTo2L2Nu_mllg->Sumw2();
  hbkg_WZTo3LNu_mllg->Sumw2();
  hbkg_WZTo2L2Q_mllg->Sumw2();
  hbkg_ZZTo2L2Nu_mllg->Sumw2();
  hbkg_ZZTo2L2Q_mllg->Sumw2();
  hbkg_ZZTo4L_mllg->Sumw2();

  hbkg_ttg_ptllg->Sumw2();
  hbkg_WWToLNuQQ_ptllg->Sumw2();
  hbkg_WWTo2L2Nu_ptllg->Sumw2();
  hbkg_WZTo3LNu_ptllg->Sumw2();
  hbkg_WZTo2L2Q_ptllg->Sumw2();
  hbkg_ZZTo2L2Nu_ptllg->Sumw2();
  hbkg_ZZTo2L2Q_ptllg->Sumw2();
  hbkg_ZZTo4L_ptllg->Sumw2();


  TCut cut = "z_mass>50 && z_charge==0 && lept0_pt>25 && lept1_pt>20 && fabs(lept0_eta)<2.4 && fabs(lept1_eta)<2.4 && gamma_pt>20";

  if (isEle) cut += "trig_Ele23_Ele12 ==1 && leptType==11";
  else cut += "trig_Mu17_Mu8 ==1 && leptType==13 && pair_dPhi>70";
  if (barrel) cut += "isEB && gamma_ChIso<2.";
  else  cut += "isEE && gamma_ChIso<1.5";
  //cut += "nselJet==0";

  TCut weight = "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF";

  
  ttg->Draw("gamma_pt >> hbkg_ttg_pt", cut*weight, "goff");
  tWWToLNuQQ->Draw("gamma_pt >> hbkg_WWToLNuQQ_pt", cut*weight, "goff");  
  tWWTo2L2Nu->Draw("gamma_pt >> hbkg_WWTo2L2Nu_pt", cut*weight, "goff");
  tWZTo3LNu->Draw("gamma_pt >> hbkg_WZTo3LNu_pt", cut*weight, "goff");
  tWZTo2L2Q->Draw("gamma_pt >> hbkg_WZTo2L2Q_pt", cut*weight, "goff");
  tZZTo2L2Nu->Draw("gamma_pt >> hbkg_ZZTo2L2Nu_pt", cut*weight, "goff");
  tZZTo2L2Q->Draw("gamma_pt >> hbkg_ZZTo2L2Q_pt", cut*weight, "goff");
  tZZTo4L->Draw("gamma_pt >> hbkg_ZZTo4L_pt", cut*weight, "goff");

  ttg->Draw("boss_mass >> hbkg_ttg_mllg", cut*weight, "goff");
  tWWToLNuQQ->Draw("boss_mass >> hbkg_WWToLNuQQ_mllg", cut*weight, "goff");
  tWWTo2L2Nu->Draw("boss_mass >> hbkg_WWTo2L2Nu_mllg", cut*weight, "goff");
  tWZTo3LNu->Draw("boss_mass >> hbkg_WZTo3LNu_mllg", cut*weight, "goff");
  tWZTo2L2Q->Draw("boss_mass >> hbkg_WZTo2L2Q_mllg", cut*weight, "goff");
  tZZTo2L2Nu->Draw("boss_mass >> hbkg_ZZTo2L2Nu_mllg", cut*weight, "goff");
  tZZTo2L2Q->Draw("boss_mass >> hbkg_ZZTo2L2Q_mllg", cut*weight, "goff");
  tZZTo4L->Draw("boss_mass >> hbkg_ZZTo4L_mllg", cut*weight, "goff");

  ttg->Draw("boss_pt >> hbkg_ttg_ptllg", cut*weight, "goff");
  tWWToLNuQQ->Draw("boss_pt >> hbkg_WWToLNuQQ_ptllg", cut*weight, "goff");
  tWWTo2L2Nu->Draw("boss_pt >> hbkg_WWTo2L2Nu_ptllg", cut*weight, "goff");
  tWZTo3LNu->Draw("boss_pt >> hbkg_WZTo3LNu_ptllg", cut*weight, "goff");
  tWZTo2L2Q->Draw("boss_pt >> hbkg_WZTo2L2Q_ptllg", cut*weight, "goff");
  tZZTo2L2Nu->Draw("boss_pt >> hbkg_ZZTo2L2Nu_ptllg", cut*weight, "goff");
  tZZTo2L2Q->Draw("boss_pt >> hbkg_ZZTo2L2Q_ptllg", cut*weight, "goff");
  tZZTo4L->Draw("boss_pt >> hbkg_ZZTo4L_ptllg", cut*weight, "goff");

  float xs_ttg = 0.635E+03; //err = 2.087e-01
  float xs_WWToLNuQQ = 50.E+03; //err = 3.475e+01
  float xs_WWTo2L2Nu = 12.178E+03; //err = 0, not known
  float xs_WZTo3LNu = 4.43E+03; //err = 7.419
  float xs_WZTo2L2Q = 5.631E+03; //err = 2.071e+01
  float xs_ZZTo2L2Nu = 0.564; //err = 0, not know
  float xs_ZZTo2L2Q = 3.25E+03; // err = 9.914
  float xs_ZZTo4L = 1.256E+03 ; //err = 2.271

  TH1F *htotwei_ttg = (TH1F*) fttg->Get("hntotweight");
  TH1F *htotwei_WWToLNuQQ = (TH1F*) fWWToLNuQQ->Get("hntotweight");
  TH1F *htotwei_WWTo2L2Nu = (TH1F*) fWWTo2L2Nu->Get("hntotweight");
  TH1F *htotwei_WZTo3LNu = (TH1F*) fWZTo3LNu->Get("hntotweight");
  TH1F *htotwei_WZTo2L2Q = (TH1F*) fWZTo2L2Q->Get("hntotweight");
  TH1F *htotwei_ZZTo2L2Nu = (TH1F*) fZZTo2L2Nu->Get("hntotweight");
  TH1F *htotwei_ZZTo2L2Q = (TH1F*) fZZTo2L2Q->Get("hntotweight");
  TH1F *htotwei_ZZTo4L = (TH1F*) fZZTo4L->Get("hntotweight");

  float lumi = 35.9;
  float scale_ttg = lumi * xs_ttg/htotwei_ttg->Integral();
  float scale_wwtolnuqq = lumi * xs_WWToLNuQQ/htotwei_WWToLNuQQ->Integral();
  float scale_wwto2lnu = lumi *xs_WWTo2L2Nu/htotwei_WWTo2L2Nu->Integral();
  float scale_wzto3lnu = lumi *xs_WZTo3LNu/htotwei_WZTo3LNu->Integral();
  float scale_wzto2l2q = lumi *xs_WZTo2L2Q/htotwei_WZTo2L2Q->Integral();
  float scale_zzto2l2nu = lumi *xs_ZZTo2L2Nu/htotwei_ZZTo2L2Nu->Integral();
  float scale_zzto2l2q = lumi *xs_ZZTo2L2Q/htotwei_ZZTo2L2Q->Integral();
  float scale_zzto4l = lumi *xs_ZZTo4L/htotwei_ZZTo4L->Integral();

  hbkg_ttg_pt->Scale(scale_ttg);
  hbkg_WWToLNuQQ_pt->Scale(scale_wwtolnuqq);
  hbkg_WWTo2L2Nu_pt->Scale(scale_wwto2lnu);
  hbkg_WZTo3LNu_pt->Scale(scale_wzto3lnu);
  hbkg_WZTo2L2Q_pt->Scale(scale_wzto2l2q);
  hbkg_ZZTo2L2Nu_pt->Scale(scale_zzto2l2nu);
  hbkg_ZZTo2L2Q_pt->Scale(scale_zzto2l2q);
  hbkg_ZZTo4L_pt->Scale(scale_zzto4l);

  hbkg_ttg_mllg->Scale(scale_ttg);
  hbkg_WWToLNuQQ_mllg->Scale(scale_wwtolnuqq);
  hbkg_WWTo2L2Nu_mllg->Scale(scale_wwto2lnu);
  hbkg_WZTo3LNu_mllg->Scale(scale_wzto3lnu);
  hbkg_WZTo2L2Q_mllg->Scale(scale_wzto2l2q);
  hbkg_ZZTo2L2Nu_mllg->Scale(scale_zzto2l2nu);
  hbkg_ZZTo2L2Q_mllg->Scale(scale_zzto2l2q);
  hbkg_ZZTo4L_mllg->Scale(scale_zzto4l);

  hbkg_ttg_ptllg->Scale(scale_ttg);
  hbkg_WWToLNuQQ_ptllg->Scale(scale_wwtolnuqq);
  hbkg_WWTo2L2Nu_ptllg->Scale(scale_wwto2lnu);
  hbkg_WZTo3LNu_ptllg->Scale(scale_wzto3lnu);
  hbkg_WZTo2L2Q_ptllg->Scale(scale_wzto2l2q);
  hbkg_ZZTo2L2Nu_ptllg->Scale(scale_zzto2l2nu);
  hbkg_ZZTo2L2Q_ptllg->Scale(scale_zzto2l2q);
  hbkg_ZZTo4L_ptllg->Scale(scale_zzto4l);

  hbkg_ttg_ptllg->Print("all");
  hbkg_ttg_pt->Add(hbkg_WWToLNuQQ_pt);
  hbkg_ttg_pt->Add(hbkg_WWTo2L2Nu_pt);
  hbkg_ttg_pt->Add(hbkg_WZTo3LNu_pt);
  hbkg_ttg_pt->Add(hbkg_WZTo2L2Q_pt);
  hbkg_ttg_pt->Add(hbkg_ZZTo2L2Nu_pt);
  hbkg_ttg_pt->Add(hbkg_ZZTo2L2Q_pt);
  hbkg_ttg_pt->Add(hbkg_ZZTo4L_pt);

  hbkg_ttg_mllg->Add(hbkg_WWToLNuQQ_mllg);
  hbkg_ttg_mllg->Add(hbkg_WWTo2L2Nu_mllg);
  hbkg_ttg_mllg->Add(hbkg_WZTo3LNu_mllg);
  hbkg_ttg_mllg->Add(hbkg_WZTo2L2Q_mllg);
  hbkg_ttg_mllg->Add(hbkg_ZZTo2L2Nu_mllg);
  hbkg_ttg_mllg->Add(hbkg_ZZTo2L2Q_mllg);
  hbkg_ttg_mllg->Add(hbkg_ZZTo4L_mllg);

  hbkg_ttg_ptllg->Add(hbkg_WWToLNuQQ_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_WWTo2L2Nu_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_WZTo3LNu_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_WZTo2L2Q_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_ZZTo2L2Nu_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_ZZTo2L2Q_ptllg);
  hbkg_ttg_ptllg->Add(hbkg_ZZTo4L_ptllg);

  hbkg_ttg_ptllg->Print("all");

  TString outname = "TTG_VV_bkg";
  if (barrel) outname += "_EBpho";
  else outname += "_EEpho";
  if (isEle) outname += "_eleChan";
  else outname += "_muChan";
  outname += ".root";

  TFile *fout = new TFile(outname, "recreate");
  fout->cd();

  hbkg_ttg_pt->Write("hbkg_ttg_vv_ptg");
  hbkg_ttg_mllg->Write("hbkg_ttg_vv_mllg");
  hbkg_ttg_ptllg->Write("hbkg_ttg_vv_ptllg");

  fout->Write();
  fout->Close();
  
}
