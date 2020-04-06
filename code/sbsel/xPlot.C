#include "TFile.h"
#include "TH1.h"

void xPlot() {

  TFile *fqcd = new TFile("input/bkgTemplate_summer16_420_BDT_Up6000_5VarCorr_job_summer16_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8.root", "read");
  //TFile *fqcd = new TFile("input/bkgTemplate_summer16_420_BDT_Up6000_5VarCorr_job_summer16_GJets_Pt_15To6000_TuneCUETP8M1_pythia8.root", "read");
  //TFile *fdy = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_etawei_isr.root", "read");
  TFile *fdy = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_matching_dr0p2.root", "read");
  TFile *fdy_sb = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_etawei.root", "read");
  TFile *fda = new TFile("../ana/minitrees/DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_etawei.root", "read");
  
  TTree *tqcd = (TTree*) fqcd->Get("outtree");
  TTree *tdy = (TTree*) fdy->Get("outtree");
  TTree *tdy_sb = (TTree*) fdy_sb->Get("outtree");
  TTree *tda = (TTree*) fda->Get("outtree");

  TH1F *hsr_dy = new TH1F("hsr_dy", "", 10, -1, 1);
  TH1F *hsb_dy = new TH1F("hsb_dy", "", 10, -1, 1);
  TH1F *hsb_qcd = new TH1F("hsb_qcd", "", 10, -1, 1);
  TH1F *hsb_da = new TH1F("hsb_da", "", 10, -1, 1);

  TH1F *hsb_dy_mu = new TH1F("hsb_dy_mu", "", 10, -1, 1);

  hsb_da->Sumw2();
  
  tdy->Draw("gamma_ssmva >> hsr_dy", "(isEE==1 && gamma_ChIso<1.5 && gamma_pt>15 && gamma_pt<20 && (leptType==11 || leptType==13))*puweigj_65nb*genWeight", "goff");
  tdy->Draw("gamma_ssmva >> hsb_dy", "(isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso<13. && ((leptType==11 && trig_Ele23_Ele12==1) || (leptType==13 && trig_Mu17_Mu8==1)))*puweigj_65nb*genWeight", "goff");
  tda->Draw("gamma_ssmva >> hsb_da", "(isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso<13. && ((leptType==11 && trig_Ele23_Ele12==1)||(leptType==13 && trig_Mu17_Mu8==1)))", "goff");
  tqcd->Draw("gamma_ssmva >> hsb_qcd", "(isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso<13.)*puwei*genwei", "goff");

  tdy->Draw("gamma_ssmva >> hsb_dy_mu", "(isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso<13. && leptType==13 && trig_Mu17_Mu8==1)*puweigj_65nb*genWeight", "goff");

  hsb_qcd->Print();

  hsr_dy_err = (TH1F*) hsr_dy->Clone();
  hsb_dy_err = (TH1F*) hsb_dy->Clone();
  hsb_qcd_err = (TH1F*) hsb_qcd->Clone();
  hsb_dy_mu_err = (TH1F*) hsb_dy_mu->Clone();
  
  hsr_dy_err->Sumw2();
  hsb_dy_err->Sumw2();
  hsb_qcd_err->Sumw2();
  hsb_dy_mu_err->Sumw2();
  

  hsr_dy->SetLineColor(2);
  hsb_dy->SetLineColor(4);
  hsb_qcd->SetLineColor(kGreen);

  hsr_dy->SetLineWidth(2);
  hsb_dy->SetLineWidth(2);
  hsb_qcd->SetLineWidth(2);

  hsr_dy_err->SetLineColor(2);
  hsr_dy_err->SetLineWidth(2);
  hsb_dy_err->SetLineColor(4);
  hsb_dy_err->SetLineWidth(2);
  hsb_qcd_err->SetLineColor(kGreen);
  hsb_qcd_err->SetLineWidth(2);

  hsr_dy_err->SetMarkerSize(0);
  hsb_dy_err->SetMarkerSize(0);
  hsb_qcd_err->SetMarkerSize(0);

  hsb_dy_mu->SetLineColor(kViolet);
  hsb_dy_mu->SetLineWidth(2);
  hsb_dy_mu_err->SetLineColor(kViolet);
  hsb_dy_mu_err->SetLineWidth(2);
  hsb_dy_mu_err->SetMarkerSize(0);
  

  hsr_dy_err->Scale(1./hsr_dy_err->Integral());
  hsb_dy_err->Scale(1./hsb_dy_err->Integral());
  hsb_qcd_err->Scale(1./hsb_qcd_err->Integral());
  hsb_dy_mu_err->Scale(1./hsb_dy_mu_err->Integral());



  hsr_dy->Scale(1./hsr_dy->Integral());
  hsb_dy->Scale(1./hsb_dy->Integral());
  hsb_da->Scale(1./hsb_da->Integral());
  hsb_qcd->Scale(1./hsb_qcd->Integral());
  hsb_dy_mu->Scale(1./hsb_dy_mu->Integral());


  gStyle->SetOptStat(0);
  
  c = new TCanvas();
  c->cd();
  c->SetLogy();
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  //pad1->SetBottomMargin(0.02);
  //pad1->Draw();
  //pad1->cd();
  //pad1->SetLogy();
  hsr_dy->GetYaxis()->SetTitle("a.u.");
  hsr_dy->GetXaxis()->SetTitle("photon MVA");
  hsr_dy->SetMaximum(0.7);
  hsr_dy->SetMinimum(1e-03);
  hsr_dy->Draw("hist");
  hsb_dy->Draw("histsame");
  hsb_da->Draw("esame");
  //hsb_qcd->Draw("histsame");
  hsr_dy_err->Draw("esames");
  hsb_dy_err->Draw("eSames");
  //hsb_qcd_err->Draw("eSames");

  //hsr_dy->GetXaxis()->SetLabelSize(0);
  TLegend *lg = new TLegend(0.28, 0.25, 0.65, 0.5);
  lg->AddEntry(hsr_dy, "signal region (DYJet)", "f");
  lg->AddEntry(hsb_dy, "7 < charged iso < 13 (SB) (DYJet)", "f");
  lg->AddEntry(hsb_da, "7 < charged iso < 13 (SB) (data)", "pe");
  //lg->AddEntry(hsb_qcd, "7 < charged iso < 13 (SB) (QCD)", "f");
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  lg->Draw();
  
  /*
  c->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();

  TH1F *hratio = (TH1F*) hsr_dy->Clone();
  hratio->Divide(hsb_dy);

  hratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hratio->GetYaxis()->SetTitle("SR/SB");
  hratio->GetXaxis()->SetTitle("photon MVA");
  hratio->SetTitleSize(0.14, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0.8, 1.4);
  hratio->GetYaxis()->SetNdivisions(8);

  hratio->Draw("");
  //hratio->Print("all");
  */

  /*
  c1 = new TCanvas();
  c1->cd();
  c1->SetLogy();
  //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  //pad1->SetBottomMargin(0.02);
  //pad1->Draw();
  //pad1->cd();
  //pad1->SetLogy();
  hsb_dy->GetYaxis()->SetTitle("a.u.");
  hsb_dy->GetXaxis()->SetTitle("photon MVA");
  hsb_dy->SetMaximum(0.7);
  hsb_dy->SetMinimum(1e-03);
  hsb_dy->Draw("hist");
  hsb_dy_mu->Draw("histsame");
  hsb_dy_err->Draw("esames");
  hsb_dy_mu_err->Draw("esames");
  */
  
  //c->SaveAs("plots/SR_SB_Pt15to20_ele_EE_data_DYJet_RemoveEle_dR0p2.pdf");



}
