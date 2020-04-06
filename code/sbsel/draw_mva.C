#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

void draw_mva(bool barrel = true) {

  TFile *fmc = TFile::Open("../ana/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr.root");

  TTree *tmc = (TTree*) fmc->Get("outtree");

  TH1F *hmva_pt15to20 = new TH1F("hmva_pt15to20", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt20to25 = new TH1F("hmva_pt20to25", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt25to35 = new TH1F("hmva_pt25to35", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt35to45 = new TH1F("hmva_pt35to45", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt45to65 = new TH1F("hmva_pt45to65", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt65to85 = new TH1F("hmva_pt65to85", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt85to120 = new TH1F("hmva_pt85to120", "hmva_pt15to20", 10, -1, 1);
  TH1F *hmva_pt120to1000 = new TH1F("hmva_pt120to1000", "hmva_pt15to20", 10, -1, 1);

  TH1F *hmva_ele = new TH1F("hmva_ele", "", 10, -1, 1);
  TH1F *hmva_mu = new TH1F("hmva_mu", "", 10, -1, 1);

  TCut cut_ele = "z_charge==0 && leptType==11 && trig_Ele23_Ele12==1 && lept0_pt>25 && lept1_pt>20 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4";
  TCut cut_mu = "z_charge==0 && leptType==13 && trig_Mu17_Mu8==1 && lept0_pt>25 && lept1_pt>20 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4";

  if (barrel) {
    cut_ele += "isEB && gamma_ChIso<2";
    cut_mu += "isEB && gamma_ChIso<2";
  }
  else {
    cut_ele += "isEE && gamma_ChIso<1.5";
    cut_mu += "isEE && gamma_ChIso<1.5";
  }

  tmc->Draw("gamma_ssmva >> hmva_ele", cut_ele);
  tmc->Draw("gamma_ssmva >> hmva_mu", cut_mu);

  //cout << hmva_ele->Integral() << endl;

  tmc->Draw("gamma_ssmva >> hmva_pt15to20", cut_ele + "gamma_pt>15 && gamma_pt<20");
  tmc->Draw("gamma_ssmva >> hmva_pt20to25", cut_ele + "gamma_pt>20 && gamma_pt<25");
  tmc->Draw("gamma_ssmva >> hmva_pt25to35", cut_ele + "gamma_pt>25 && gamma_pt<35");
  tmc->Draw("gamma_ssmva >> hmva_pt35to45", cut_ele + "gamma_pt>35 && gamma_pt<45");
  tmc->Draw("gamma_ssmva >> hmva_pt45to65", cut_ele + "gamma_pt>45 && gamma_pt<65");
  tmc->Draw("gamma_ssmva >> hmva_pt65to85", cut_ele + "gamma_pt>65 && gamma_pt<85");
  tmc->Draw("gamma_ssmva >> hmva_pt85to120", cut_ele + "gamma_pt>85 && gamma_pt<1000");
  tmc->Draw("gamma_ssmva >> hmva_pt120to1000", cut_ele + "gamma_pt>120 && gamma_pt<1000");

  hmva_ele->Scale(1./hmva_ele->Integral());
  hmva_mu->Scale(1./hmva_mu->Integral());

  hmva_pt15to20->Scale(1./hmva_pt15to20->Integral());
  hmva_pt20to25->Scale(1./hmva_pt20to25->Integral());
  hmva_pt25to35->Scale(1./hmva_pt25to35->Integral());
  hmva_pt35to45->Scale(1./hmva_pt35to45->Integral());
  hmva_pt45to65->Scale(1./hmva_pt45to65->Integral());
  hmva_pt65to85->Scale(1./hmva_pt65to85->Integral());
  hmva_pt85to120->Scale(1./hmva_pt85to120->Integral());
  hmva_pt120to1000->Scale(1./hmva_pt120to1000->Integral());

  hmva_ele->SetLineColor(2);
  hmva_mu->SetLineColor(4);

  hmva_pt15to20->SetLineColor(kRed);
  hmva_pt20to25->SetLineColor(kRed-6);
  hmva_pt25to35->SetLineColor(kOrange);
  hmva_pt35to45->SetLineColor(kBlue);
  hmva_pt45to65->SetLineColor(kGreen+1);
  hmva_pt65to85->SetLineColor(kBlue-8);
  hmva_pt85to120->SetLineColor(kCyan-2);
  hmva_pt120to1000->SetLineColor(kViolet);

  hmva_pt15to20->SetLineWidth(2);
  hmva_pt20to25->SetLineWidth(2);
  hmva_pt25to35->SetLineWidth(2);
  hmva_pt35to45->SetLineWidth(2);
  hmva_pt45to65->SetLineWidth(2);
  hmva_pt65to85->SetLineWidth(2);
  hmva_pt85to120->SetLineWidth(2);
  hmva_pt120to1000->SetLineWidth(2);

  gStyle->SetOptStat(0);
  
  c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  hmva_ele->GetYaxis()->SetTitle("a.u");
  hmva_ele->GetXaxis()->SetTitle("photon MVA");
  hmva_ele->Draw();
  hmva_mu->Draw("same");
  TLegend *lg = new TLegend(0.3, 0.4, 0.48, 0.53);
  lg->AddEntry(hmva_ele, "ele", "f");
  lg->AddEntry(hmva_mu, "muon", "f");
  lg->SetBorderSize(0);
  lg->Draw();

  c2 = new TCanvas("c2", "c2", 650, 650);
  c2->cd();
  hmva_pt15to20->GetYaxis()->SetTitle("a.u");
  hmva_pt15to20->GetXaxis()->SetTitle("photon MVA");
  hmva_pt15to20->GetYaxis()->SetRangeUser(0., 0.4);
  hmva_pt15to20->Draw();
  hmva_pt20to25->Draw("same");
  hmva_pt25to35->Draw("same");
  hmva_pt35to45->Draw("same");
  hmva_pt45to65->Draw("same");
  hmva_pt65to85->Draw("same");
  hmva_pt85to120->Draw("same");
  hmva_pt120to1000->Draw("same");
  TLegend *lg1 = new TLegend(0.25, 0.44, 0.5, 0.83);
  lg1->AddEntry(hmva_pt15to20, "15 < p_{T}^{#gamma} < 20", "f");
  lg1->AddEntry(hmva_pt20to25, "20 < p_{T}^{#gamma} < 25", "f");
  lg1->AddEntry(hmva_pt25to35, "25 < p_{T}^{#gamma} < 35", "f");
  lg1->AddEntry(hmva_pt35to45, "35 < p_{T}^{#gamma} < 45", "f");
  lg1->AddEntry(hmva_pt45to65, "45 < p_{T}^{#gamma} < 65", "f");
  lg1->AddEntry(hmva_pt65to85, "65 < p_{T}^{#gamma} < 85", "f");
  lg1->AddEntry(hmva_pt85to120, "85 < p_{T}^{#gamma} < 120", "f");
  lg1->AddEntry(hmva_pt120to1000, "120 < p_{T}^{#gamma} < 1000", "f");
  lg1->SetBorderSize(0);
  lg1->SetTextFont(42);
  lg1->Draw();

  
  

}
