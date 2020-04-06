#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

#include <iostream>

using namespace std;

void getBkg_njet(bool isEle = 1) {

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

  int njet = 7;

  TH1F *hbkg_ttg_vv_njet_ele = new TH1F("hbkg_ttg_vv_njet_ele", "hbkg_ttg_v of njet", njet, 0, njet);
  TH1F *hbkg_ttg_vv_njet_mu = new TH1F("hbkg_ttg_vv_njet_mu", "hbkg_ttg_v of njet", njet, 0, njet);

  hbkg_ttg_vv_njet_ele->Sumw2();
  hbkg_ttg_vv_njet_mu->Sumw2();

  TCut cut = "z_mass>50 && z_charge==0 && lept0_pt>25 && lept1_pt>20 && fabs(lept0_eta)<2.4 && fabs(lept1_eta)<2.4 && gamma_pt>20";
  //if (isEle) cut += "trig_Ele23_Ele12 ==1 && leptType==11";
  //else cut += "trig_Mu17_Mu8 ==1 && leptType==13 && pair_dPhi>70";
  cut += "(isEB && gamma_ChIso<2.) || (isEE && gamma_ChIso<1.5)";

  TCut cutele = "trig_Ele23_Ele12 ==1 && leptType==11";
  TCut cutmu = "trig_Mu17_Mu8 ==1 && leptType==13 && pair_dPhi>70";

  TCut weight = "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF";

  for (int i = 0; i < njet; i++) {

    TH1F *hpt_ttg_ele = new TH1F("hpt_ttg_ele", "pt", 1, 15, 1000);
    TH1F *hpt_wwlnuqq_ele = new TH1F("hpt_wwlnuqq_ele", "pt", 1, 15, 1000);
    TH1F *hpt_ww2l2nu_ele = new TH1F("hpt_ww2l2nu_ele", "pt", 1, 15, 1000);
    TH1F *hpt_wz3lnu_ele  = new TH1F("hpt_wz3lnu_ele", "pt", 1, 15, 1000);
    TH1F *hpt_wz2l2q_ele  = new TH1F("hpt_wz2l2q_ele", "pt", 1, 15, 1000);
    TH1F *hpt_zz2l2nu_ele = new TH1F("hpt_zz2l2nu_ele", "pt", 1, 15, 1000);
    TH1F *hpt_zz2l2q_ele  = new TH1F("hpt_zz2l2q_ele", "pt", 1, 15, 1000);
    TH1F *hpt_zz4l_ele    = new TH1F("hpt_zz4l_ele", "pt", 1, 15, 1000);

    TH1F *hpt_ttg_mu = new TH1F("hpt_ttg_mu", "pt", 1, 15, 1000);
    TH1F *hpt_wwlnuqq_mu = new TH1F("hpt_wwlnuqq_mu", "pt", 1, 15, 1000);
    TH1F *hpt_ww2l2nu_mu = new TH1F("hpt_ww2l2nu_mu", "pt", 1, 15, 1000);
    TH1F *hpt_wz3lnu_mu  = new TH1F("hpt_wz3lnu_mu", "pt", 1, 15, 1000);
    TH1F *hpt_wz2l2q_mu  = new TH1F("hpt_wz2l2q_mu", "pt", 1, 15, 1000);
    TH1F *hpt_zz2l2nu_mu = new TH1F("hpt_zz2l2nu_mu", "pt", 1, 15, 1000);
    TH1F *hpt_zz2l2q_mu  = new TH1F("hpt_zz2l2q_mu", "pt", 1, 15, 1000);
    TH1F *hpt_zz4l_mu    = new TH1F("hpt_zz4l_mu", "pt", 1, 15, 1000);

    TCut cutjet;
    if (i==0) cutjet= "nselJet==0";
    if (i==1) cutjet= "nselJet==1 && jetPt[0]>30. && jetPt[0]<50.";
    if (i==2) cutjet= "nselJet==1 && jetPt[0]>=50.";
    if (i==3) cutjet= "nselJet==2 && jetPt[1]>30. && jetPt[1]<50.";
    if (i==4) cutjet= "nselJet==2 && jetPt[1]>=50.";
    if (i==5) cutjet= "nselJet>=3 && jetPt[2]>30. && jetPt[2]<50.";
    if (i==6) cutjet= "nselJet>=3 && jetPt[2]>=50.";


    ttg->Draw("gamma_pt >> hpt_ttg_ele", (cut+cutele+cutjet)*weight, "goff");
    tWWToLNuQQ->Draw("gamma_pt >> hpt_wwlnuqq_ele", (cut+cutele+cutjet)*weight, "goff");
    tWWTo2L2Nu->Draw("gamma_pt >> hpt_ww2l2nu_ele", (cut+cutele+cutjet)*weight, "goff");
    tWZTo3LNu->Draw("gamma_pt >> hpt_wz3lnu_ele", (cut+cutele+cutjet)*weight, "goff");
    tWZTo2L2Q->Draw("gamma_pt >> hpt_wz2l2q_ele", (cut+cutele+cutjet)*weight, "goff");
    tZZTo2L2Nu->Draw("gamma_pt >> hpt_zz2l2nu_ele", (cut+cutele+cutjet)*weight, "goff");
    tZZTo2L2Q->Draw("gamma_pt >> hpt_zz2l2q_ele", (cut+cutele+cutjet)*weight, "goff");
    tZZTo4L->Draw("gamma_pt >> hpt_zz4l_ele", (cut+cutele+cutjet)*weight, "goff");

    ttg->Draw("gamma_pt >> hpt_ttg_mu", (cut+cutmu+cutjet)*weight, "goff");
    tWWToLNuQQ->Draw("gamma_pt >> hpt_wwlnuqq_mu", (cut+cutmu+cutjet)*weight, "goff");
    tWWTo2L2Nu->Draw("gamma_pt >> hpt_ww2l2nu_mu", (cut+cutmu+cutjet)*weight, "goff");
    tWZTo3LNu->Draw("gamma_pt >> hpt_wz3lnu_mu", (cut+cutmu+cutjet)*weight, "goff");
    tWZTo2L2Q->Draw("gamma_pt >> hpt_wz2l2q_mu", (cut+cutmu+cutjet)*weight, "goff");
    tZZTo2L2Nu->Draw("gamma_pt >> hpt_zz2l2nu_mu", (cut+cutmu+cutjet)*weight, "goff");
    tZZTo2L2Q->Draw("gamma_pt >> hpt_zz2l2q_mu", (cut+cutmu+cutjet)*weight, "goff");
    tZZTo4L->Draw("gamma_pt >> hpt_zz4l_mu", (cut+cutmu+cutjet)*weight, "goff");

    cut.Print();
    cutjet.Print();

    float xs_ttg = 0.635E+03; //err = 2.087e-01
    float xs_WWToLNuQQ = 50.E+03; //err = 3.475e+01
    float xs_WWTo2L2Nu = 12.178E+03; //err = 0, not known
    float xs_WZTo3LNu = 4.43E+03; //err = 7.419
    float xs_WZTo2L2Q = 5.631E+03; //err = 2.071e+01
    float xs_ZZTo2L2Nu = 0.564; //err = 0, not know
    float xs_ZZTo2L2Q = 3.25E+03 ; // err = 9.914
    float xs_ZZTo4L = 1.256E+03  ; //err = 2.271
    
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
    
    hpt_ttg_ele->Scale(scale_ttg);
    hpt_wwlnuqq_ele->Scale(scale_wwtolnuqq);
    hpt_ww2l2nu_ele->Scale(scale_wwto2lnu);
    hpt_wz3lnu_ele->Scale(scale_wzto3lnu);
    hpt_wz2l2q_ele->Scale(scale_wzto2l2q);
    hpt_zz2l2nu_ele->Scale(scale_zzto2l2nu);
    hpt_zz2l2q_ele->Scale(scale_zzto2l2q);
    hpt_zz4l_ele->Scale(scale_zzto4l);

    hpt_ttg_mu->Scale(scale_ttg);
    hpt_wwlnuqq_mu->Scale(scale_wwtolnuqq);
    hpt_ww2l2nu_mu->Scale(scale_wwto2lnu);
    hpt_wz3lnu_mu->Scale(scale_wzto3lnu);
    hpt_wz2l2q_mu->Scale(scale_wzto2l2q);
    hpt_zz2l2nu_mu->Scale(scale_zzto2l2nu);
    hpt_zz2l2q_mu->Scale(scale_zzto2l2q);
    hpt_zz4l_mu->Scale(scale_zzto4l);

    //bkg_ttg[i] = hpt_ttg->Integral();
    //bkg_WWToLNuQQ[i] = hpt_wwlnuqq->Integral();
    //bkg_WWTo2L2Nu[i] = hpt_ww2l2nu->Integral();
    //bkg_WZTo3LNu[i] = hpt_wz3lnu->Integral();
    //bkg_WZTo2L2Q[i] = hpt_wz2l2q->Integral();
    //bkg_ZZTo2L2Nu[i] = hpt_zz2l2nu->Integral();
    //bkg_ZZTo2L2Q[i] = hpt_zz2l2q->Integral();
    //bkg_ZZTo4L[i] = hpt_zz4l->Integral();

    float totbkg_ele = hpt_ttg_ele->Integral() + hpt_wwlnuqq_ele->Integral() + hpt_ww2l2nu_ele->Integral() 
                   + hpt_wz3lnu_ele->Integral() + hpt_wz2l2q_ele->Integral() 
                   + hpt_zz2l2nu_ele->Integral() + hpt_zz2l2q_ele->Integral() + hpt_zz4l_ele->Integral();

    float toterr_ele = pow(hpt_ttg_ele->GetBinError(1),2) + pow(hpt_wwlnuqq_ele->GetBinError(1),2) + pow(hpt_ww2l2nu_ele->GetBinError(1),2)
                 + pow(hpt_wz3lnu_ele->GetBinError(1),2) +pow(hpt_wz2l2q_ele->GetBinError(1),2)
      + pow(hpt_zz2l2nu_ele->GetBinError(1),2) + pow(hpt_zz2l2q_ele->GetBinError(1),2) +pow(hpt_zz4l_ele->GetBinError(1),2);

    
    float totbkg_mu = hpt_ttg_mu->Integral() + hpt_wwlnuqq_mu->Integral() + hpt_ww2l2nu_mu->Integral() 
                   + hpt_wz3lnu_mu->Integral() + hpt_wz2l2q_mu->Integral() 
                   + hpt_zz2l2nu_mu->Integral() + hpt_zz2l2q_mu->Integral() + hpt_zz4l_mu->Integral();

    float toterr_mu = pow(hpt_ttg_mu->GetBinError(1),2) + pow(hpt_wwlnuqq_mu->GetBinError(1),2) + pow(hpt_ww2l2nu_mu->GetBinError(1),2)
                 + pow(hpt_wz3lnu_mu->GetBinError(1),2) +pow(hpt_wz2l2q_mu->GetBinError(1),2)
      + pow(hpt_zz2l2nu_mu->GetBinError(1),2) + pow(hpt_zz2l2q_mu->GetBinError(1),2) +pow(hpt_zz4l_mu->GetBinError(1),2);

    hbkg_ttg_vv_njet_ele->SetBinContent(i+1, totbkg_ele);
    hbkg_ttg_vv_njet_ele->SetBinError(i+1, sqrt(toterr_ele));

    hbkg_ttg_vv_njet_mu->SetBinContent(i+1, totbkg_mu);
    hbkg_ttg_vv_njet_mu->SetBinError(i+1, sqrt(toterr_mu));

    hpt_ttg_ele->Delete();
    hpt_wwlnuqq_ele->Delete();
    hpt_ww2l2nu_ele->Delete();
    hpt_wz3lnu_ele->Delete();
    hpt_wz2l2q_ele->Delete();
    hpt_zz2l2nu_ele->Delete();
    hpt_zz2l2q_ele->Delete();
    hpt_zz4l_ele->Delete();

    hpt_ttg_mu->Delete();
    hpt_wwlnuqq_mu->Delete();
    hpt_ww2l2nu_mu->Delete();
    hpt_wz3lnu_mu->Delete();
    hpt_wz2l2q_mu->Delete();
    hpt_zz2l2nu_mu->Delete();
    hpt_zz2l2q_mu->Delete();
    hpt_zz4l_mu->Delete();
}

  TString outname = "TTG_VV_bkg";
  //if (isEle) outname += "_eleChan";
  //else outname += "_muChan";
  outname += "_njet.root";

  TFile *fout = new TFile(outname, "recreate");
  fout->cd();

  hbkg_ttg_vv_njet_ele->Write();
  hbkg_ttg_vv_njet_mu->Write();

  fout->Write();
  fout->Close();
  
}
