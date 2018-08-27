#include "TFile.h"
#include "TH1.h"
#include "TTree.h"


void ZmassPlot(bool ele = true) {

  gStyle->SetOptStat(0);

  TChain *tda = new TChain("tZ");
  if (ele) {
    tda->Add("../ana/minitrees/DoubleEG_Run2016B_FebReminiAOD_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016C_FebReminiAOD_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016D_FebReminiAOD_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016E_FebReminiAOD_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016F_FebReminiAOD1_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016F_FebReminiAOD2_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016G_FebReminiAOD_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016H_FebReminiAODv2_SortedPt.root");
    tda->Add("../ana/minitrees/DoubleEG_Run2016H_FebReminiAODv3_SortedPt.root");
  }

  else {
  }

  TFile *fsig = new TFile("../ana/minitrees/ZJets_aMCatNLO_SortedPt.root", "read");
  TFile *fbkg = new TFile("../ana/minitrees/TT_Powheg_SortedPt.root", "read");
  TTree *tsig = (TTree*) fsig->Get("tZ");
  TTree *tbkg = (TTree*) fbkg->Get("tZ");

  TH1F *hZm_da = new TH1F("hZm_da", "; M_{ll} (GeV); Events/2 GeV", 50, 50, 150);
  TH1F *hZm_da_ss = new TH1F("hZm_da_ss", "; M_{ll} (GeV); Events/2 GeV", 50, 50, 150);
  TH1F *hZm_sig = new TH1F("hZm_sig", "; M_{ll} (GeV); Events/2 GeV", 50, 50, 150);
  TH1F *hZm_bkg = new TH1F("hZm_bkg", "; M_{ll} (GeV); Events/2 GeV", 50, 50, 150);

  TCut cut;
  if (ele) cut = "leptType==11 && z_charge==0 && trig_Ele23_Ele12==1";
  else cut = "leptType==13 && z_charge==0 && trig_Mu17_Mu8==1";

  TCut weight;
  if (ele) weight = "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  else weight = "puweigj*genWeight*lept0_SelSF*lept1_SelSF*lept0_trigSF*lept1_trigSF*lept_dzSF";

  tda->Draw("z_mass >> hZm_da", cut);
  tda->Draw("z_mass >> hZm_da_ss", "leptType==11 && z_charge!=0 && trig_Ele23_Ele12==1");
  tsig->Draw("z_mass >> hZm_sig", cut*weight);
  tbkg->Draw("z_mass >> hZm_bkg", cut*weight);

  //hZm_da->Add(hZm_da_ss, -1);
  hZm_da->Sumw2();
  hZm_da_ss->Sumw2();
  hZm_da_ss->SetMarkerStyle(24);
  hZm_da_ss->SetMarkerColor(2);
  hZm_da_ss->SetLineColor(2);

  cout << "event from data: " << hZm_da->Integral(-1,-1) << endl;

  TH1F *htotwei_Zg = (TH1F*) fsig->Get("hntotweight");
  double totEvwei = htotwei_Zg->Integral();

  TH1F *htotwei_tt = (TH1F*) fbkg->Get("hntotweight");
  double totwei_tt = htotwei_tt->Integral();

  float lumi = 35900.;

  float w_Zg = (lumi*5943.2/totEvwei);
  float w_tt = (lumi*87.31/totwei_tt);

  hZm_sig->Scale(w_Zg);
  hZm_bkg->Scale(w_tt);

  cout << "event from sig: " << hZm_sig->Integral(-1,-1) << endl;
  cout << "event from bkg: " << hZm_bkg->Integral(-1,-1) << endl;

  hZm_sig->SetLineColor(kOrange);
  hZm_sig->SetFillColor(kOrange);

  hZm_bkg->SetLineColor(kMagenta-7);
  hZm_bkg->SetFillColor(kMagenta-7);

  TH1F *hmc = (TH1F*) hZm_sig->Clone();
  hmc->Add(hZm_bkg);

  TH1F* hRatio = (TH1F*) hZm_da->Clone();
  hRatio->Divide(hmc);
  hRatio->SetTitleSize(0.14, "XYZ");
  hRatio->SetLabelSize(0.14, "XYZ");
  hRatio->SetTitleFont(42, "XYZ");
  hRatio->GetYaxis()->SetTitle("data/ MC");
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetXaxis()->SetTitleOffset(1.1);
  hRatio->GetYaxis()->SetRangeUser(0.8, 1.2);
  hRatio->GetYaxis()->SetNdivisions(5);

  THStack *hstack = new THStack("hstack", "");
  hstack->Add(hZm_bkg);
  hstack->Add(hZm_sig);

  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  c1->SetBottomMargin(0.15);
  //c1->SetLogy();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  hZm_da->GetYaxis()->SetRangeUser(0.9, hZm_da->GetMaximum()*10);
  hZm_da->GetXaxis()->SetLabelSize(0);
  hZm_da->Draw();
  //hZm_sig->Draw("histsame");
  hstack->Draw("histsame");
  hZm_da->Draw("same");
  hZm_da_ss->Draw("same");
  pad1->RedrawAxis();
  TLegend *lg = new TLegend(0.68, 0.68, 0.85, 0.88);
  lg->AddEntry(hZm_da, "data", "p");
  lg->AddEntry(hZm_sig, "DY", "f");
  //lg->AddEntry(h[2], "ZJets", "f");
  lg->AddEntry(hZm_bkg, "T#barT", "f");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->Draw();
  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);

  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  hRatio->GetXaxis()->SetTickLength(0.05);
  hRatio->GetXaxis()->SetTitle("M_{ll} (GeV)");
  hRatio->Draw("ep");

}
