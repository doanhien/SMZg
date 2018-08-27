#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCut.h>

#include <iostream>

//#include "tdrstyle.C"
#include "CMS_lumi.C"


using namespace std;

void histStyle(TH1F *h1) {
  h1->SetTitleSize(0.14, "XYZ");
  h1->SetLabelSize(0.14, "XYZ");
  h1->SetTitleFont(42, "XYZ");
  h1->GetYaxis()->SetTitle("data/ MC");
  //h1->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  h1->GetYaxis()->SetTitleOffset(0.55);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetRangeUser(0.6, 1.4);
  h1->GetYaxis()->SetNdivisions(5);
  //h1->Draw("ep");
}

void histStyle1 (TH1F *h) {
  h->SetFillColor(kOrange+1);
  h->SetLineColor(kOrange+1);
  //h->SetFillColor(kAzure+6);
  //h->SetLineColor(kAzure+6);
}

void histStyle2 (TH1F *h) {
  h->SetFillColor(kYellow);
  h->SetLineColor(kYellow);
  h->SetFillColor(kCyan);
  h->SetLineColor(kCyan);
}

void histStyle3 (TH1F *h) {
  h->SetFillColor(kBlue-9);
  h->SetLineColor(kBlue-9);

}

void histStyle4 (TH1F *h) {
  h->SetFillColor(kBlue-4);
  h->SetLineColor(kBlue-4);

  //h->SetFillColor(kViolet-7);
  //h->SetLineColor(kViolet-7);

}


void dataMCPlot(bool ele = true, TString hname = "hlept1_Pt") {

  gStyle->SetOptStat(0);

  cout << "get started" << endl;
  TString infname;
  if (ele) infname  = "histo/Kinematics/Ele_EB_";
  else infname  = "histo/Kinematics/Mu_EB_";

  const int nFile = 11;
  TFile *fin[nFile];
  TH1F *h[nFile];

  TH1F  *hMC = new TH1F();

  if (ele) 
    fin[0] = new TFile(infname + "DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_nodoublecount.root", "read");
  else fin[0] = new TFile(infname + "DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_nodoublecount.root", "read");
  fin[1] = new TFile(infname + "ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root", "read");
  //fin[1] = new TFile(infname + "Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root", "read");
  //fin[2] = new TFile(infname + "ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_allfsr.root", "read");
  fin[2] = new TFile(infname + "ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root", "read");
  fin[3] = new TFile(infname + "TT_Powheg_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[4] = new TFile(infname + "WWTo2L2Nu_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[5] = new TFile(infname + "WWToLNuQQ_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[6] = new TFile(infname + "WZTo3LNu_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[7] = new TFile(infname + "WZTo2L2Q_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[8] = new TFile(infname + "ZZTo2L2Nu_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[9] = new TFile(infname + "ZZTo2L2Q_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  fin[10] = new TFile(infname + "ZZTo4L_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");

  TH1F *hzm_ss;
  //cout << "done opening files" << endl;
  for (int i = 0; i < nFile; i++) {
    h[i] = (TH1F*) fin[i]->Get(hname.Data());
  }

  hzm_ss = (TH1F*) fin[0]->Get("hZmass_ss");

  for (int i = 1; i < nFile; i++) {
    if (i == 1) hMC = (TH1F*) h[i]->Clone();
    else hMC->Add(h[i]);
  }

  //hMC = (TH1F*) h[2]->Clone();  
  //hMC->Add(h[3]); 
  //hMC->Add(h[4]);
  //hMC->Add(h[5]);

  h[4]->Add(h[5]);
  h[4]->Add(h[6]);
  h[4]->Add(h[7]);
  h[4]->Add(h[8]);
  h[4]->Add(h[9]);
  h[4]->Add(h[10]);

  //cout << "add hist of all mc together" << endl;
  cout << "total entries from mc: " << hMC->Integral() << endl;
  cout << "total entries from data: " << h[0]->Integral() << endl;
  cout << "total entries from Zg: " << h[1]->Integral() << endl;
  cout << "total entries from Zj: " << h[2]->Integral() << endl;
  cout << "total entries from TT: " << h[3]->Integral() << endl;
  cout << "total entries from VV: " << h[4]->Integral() << endl;

  //cout << "histogram of ratio of data and mc" << endl;

  //if (hname == "hZmass") h[0]->Add(hzm_ss, -1);
  h[0]->Sumw2();
  hzm_ss->Sumw2();
  hzm_ss->SetMarkerStyle(33);
  hzm_ss->SetMarkerColor(kRed);
  hzm_ss->SetLineColor(kRed);

  TH1F *hRatio = (TH1F*) h[0]->Clone();
  hRatio->Divide(hMC);
  histStyle(hRatio);

  histStyle1(h[1]);
  histStyle2(h[2]);
  histStyle3(h[3]);
  histStyle4(h[4]);

  THStack *hstack = new THStack("hstack", "");
  hstack->Add(h[4]);
  hstack->Add(h[3]);
  hstack->Add(h[2]);
  hstack->Add(h[1]);

  TString xname, yname, outname;
  outname = "plots/Kinematics/PhoPreSel_ChIso/data_MC_"; 

  if (ele) {
    if (hname == "hlept1_Pt") {
      xname = "leading p_{T}^{e} (GeV)";   yname = "Events"; outname += "leadEle_Pt";}
    if (hname == "hlept2_Pt") {
      xname = "trailing p_{T}^{e} (GeV)";   yname = "Events"; outname += "trailEle_Pt";}
    if (hname == "hlept1_Eta") {
      xname = "leading #eta^{e}";   yname = "Events"; outname += "leadEle_Eta"; }
    if (hname == "hlept2_Eta") {
      xname = "trailing #eta^{e}";   yname = "Events"; outname += "trailEle_Eta";  }
    if (hname == "hZmass") {
      xname = "M_{ee} (GeV)";   yname = "Events"; outname += "Mee"; }
    if (hname == "hphoEt") {
      xname = "p_{T}^{#gamma}";   yname = "Events"; outname += "Ele_PhoEt"; }
    if (hname == "hphoEta") {
      xname = "#eta^{#gamma}";   yname = "Events"; outname += "Ele_PhoEta"; }
    if (hname == "hmllg") {
      xname = "M_{ee#gamma} (GeV)";   yname = "Events"; outname += "Meeg"; }
    if (hname == "hdRPhoLep") {
      xname = "#Delta(#gamma,e)";   yname = "Events"; outname += "dRPhoEle"; }
    if (hname == "hnVtx") {
      xname = "number of vertice";   yname = "Events"; outname += "Ele_nVtx"; }
    if (hname == "hssmva") {
      xname = "showershape MVA";   yname = "Events/0.1"; outname += "Ele_ssmva"; }
  }
  else {
    if (hname == "hlept1_Pt") {
      xname = "leading p_{T}^{#mu} (GeV)";   yname = "Events"; outname += "leadMu_Pt";}
    if (hname == "hlept2_Pt") {
      xname = "trailing p_{T}^{#mu} (GeV)";   yname = "Events"; outname += "trailMu_Pt";}
    if (hname == "hlept1_Eta") {
      xname = "leading #eta^{#mu}";   yname = "Events"; outname += "leadMu_Eta"; }
    if (hname == "hlept2_Eta") {
      xname = "trailing #eta^{#mu}";   yname = "Events"; outname += "trailMu_Eta";  }
    if (hname == "hZmass") {
      xname = "M_{#mu#mu} (GeV)";   yname = "Events"; outname += "Muu"; }
    if (hname == "hphoEt") {
      xname = "p_{T}^{#gamma}";   yname = "Events"; outname += "Mu_PhoEt"; }
    if (hname == "hphoEta") {
      xname = "#eta^{#gamma}";   yname = "Events"; outname += "Mu_PhoEta"; }
    if (hname == "hmllg") {
      xname = "M_{#mu#mu#gamma} (GeV)";   yname = "Events"; outname += "Muug"; }
    if (hname == "hdRPhoLep") {
      xname = "#Delta(#gamma,#mu)";   yname = "Events"; outname += "dRPhoMu"; }
    if (hname == "hnVtx") {
      xname = "number of vertice";   yname = "Events"; outname += "Mu_nVtx"; }
    if (hname == "hssmva") {
      xname = "showershape MVA";   yname = "Events/0.1"; outname += "Mu_ssmva"; }

  }

  //cout << "outname for saving: " << outname.Data() << endl;
  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.05);
  tx.SetTextFont(52);

  //plotting for canvas
  TCanvas* c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  c1->SetBottomMargin(0.15);
  //c1->SetLogy();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  h[0]->GetYaxis()->SetRangeUser(0.9, h[0]->GetMaximum()*30);
  if (hname.Contains("Eta") || hname.Contains("Vtx") )
      h[0]->GetYaxis()->SetRangeUser(0.9, h[0]->GetMaximum()*30);
  h[0]->GetXaxis()->SetLabelSize(0);
  h[0]->GetYaxis()->SetTitle(yname);
  h[0]->Draw();
  hstack->Draw("histsame");
  h[0]->Draw("same");
  //if (hname == "hZmass") 
  //hzm_ss->Draw("same");
  pad1->RedrawAxis();
  TLegend *lg = new TLegend(0.70, 0.64, 0.85, 0.88);
  lg->AddEntry(h[0], "data", "pe");
  lg->AddEntry(h[1], "Z#gamma", "f");
  lg->AddEntry(h[2], "DYJets", "f");
  lg->AddEntry(h[3], "T#barT", "f");
  lg->AddEntry(h[4], "VV", "f");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->Draw();
  if (infname.Contains("2016B")) tx.DrawLatex(0.3, 0.8, "Run2016B");
  else if (infname.Contains("2016C")) tx.DrawLatex(0.3, 0.8, "Run2016C");
  else if (infname.Contains("2016D")) tx.DrawLatex(0.3, 0.8, "Run2016D");
  else if (infname.Contains("2016E")) tx.DrawLatex(0.3, 0.8, "Run2016E");
  else if (infname.Contains("2016F")) tx.DrawLatex(0.3, 0.8, "Run2016F");
  else if (infname.Contains("2016G")) tx.DrawLatex(0.3, 0.8, "Run2016G");
  else if (infname.Contains("2016H")) tx.DrawLatex(0.3, 0.8, "Run2016H");

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  hRatio->GetXaxis()->SetTickLength(0.05);
  hRatio->GetXaxis()->SetTitle(xname);
  hRatio->Draw("ep");
  //if (hname == "hphoEt") hRatio->Draw("epl");

  CMS_lumi(pad1, 4, 0);

  c1->SaveAs(outname + "_TightID_TMVA420_UpTo6000_EB_genmatch_dr0p1.pdf");
  //c1->SaveAs(outname + "_TightID_TMVA420_UpTo6000_5VarCorr.pdf");

}
