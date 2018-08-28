#include "TFile.h"
#include "TH1.h"

void xPlot() {


  gStyle->SetOptStat(0);

  TFile *fin = new TFile("../ana/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr_nodoublecount.root", "read");
  //TFile *fin = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_allfsr.root", "read");
  TFile *f1 = new TFile("../ana/minitrees/Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_addfsr_dr0p1.root", "read");

  TTree *t = (TTree*) fin->Get("outtreeGen");
  TTree *t1 = (TTree*) f1->Get("outtreeGen");
  TTree *tnew = (TTree*) fin->Get("outtree");
  TTree *told = (TTree*) f1->Get("outtree");


  TH1F *hgen_elept = new TH1F("hgen_elept", "", 50, 0, 150);
  TH1F *hgen_mupt = new TH1F("hgen_mupt", "", 50, 0, 150);

  TH1F *hgen_elept_nofsr = new TH1F("hgen_elept_nofsr", "", 50, 0, 150);
  TH1F *hgen_mupt_nofsr = new TH1F("hgen_mupt_nofsr", "", 50, 0, 150);

  TH1F *hgenZee = new TH1F("hgenZee", "", 70, 0, 210);
  TH1F *hgenZmm = new TH1F("hgenZmm", "", 70, 0, 210);

  TH1F *hmva_zg_old = new TH1F("hmva_zg_old", "photon mva", 20, -1, 1);
  TH1F *hmva_zg_new = new TH1F("hmva_zg_new", "photon mva", 20, -1, 1);

  TH1F *hdral = new TH1F("hdral", "hdral", 200, 0, 4);
  TH1F *hpho_fsr = new TH1F("hpho_fsr", "hpho_fs", 50, 0, 10);
  TH1F *hdrlfsr = new TH1F("hdrlfsr", "hdrlfsr", 100, 0, 5);

  //t->Draw("gendRPhoLep2 >> hdral", "genPhoEt>0. && (ngenEle==2 || ngenMu==2)");
  t->Draw("dRlhePhoEle2 >> hdral", "lheEle2>0. && lhePho1>0.");
  t->Draw("fsr_gen_pho_pt >> hpho_fsr", "fsr_gen_pho_pt>0. && no_fsr_pho>0.");
  t->Draw("fsr_pho_lep_dr1 >> hdrlfsr", "fsr_gen_pho_pt>0. && no_fsr_pho>0. && (ngenEle==2 || ngenMu==2)");
   
  t->Draw("genlepPt[1] >> hgen_elept", "ngenEle==2 && genZm>0. && genlepPt[0]>0 && genlepPt[1]>0", "goff");
  t->Draw("genlepPt[1] >> hgen_mupt", "ngenMu==2 && genZm>0. && genlepPt[0]>0 && genlepPt[1]>0", "goff");

  //hgen_elept->Scale(1./hgen_elept->Integral());
  //hgen_mupt->Scale(1./hgen_mupt->Integral());
  
  //t->Draw("lheEle1 >> hgen_elept", "", "goff");
  //t->Draw("lheMu1 >> hgen_mupt", "", "goff");

  t1->Draw("genlepPt[0] >> hgen_elept_nofsr", "ngenEle==2 && genZm>0", "goff");
  t1->Draw("genlepPt[0] >> hgen_mupt_nofsr", "ngenMu==2 && genZm>0", "goff");

  t->Draw("mcZm >> hgenZee", "ngenEle==2 && genZm>0", "goff");
  t->Draw("mcZm >> hgenZmm", "ngenMu==2 && genZm>0", "goff");

  tnew->Draw("gamma_ssmva >> hmva_zg_new", "isEE && leptType==11");
  told->Draw("gamma_ssmva >> hmva_zg_old", "isEE && leptType==11");

  
  hgen_elept->Print();
  hgen_mupt->Print();
  //hgenZee->Print();
  //hgenZmm->Print();

  
  hgen_elept->SetLineColor(2);
  hgen_elept->SetLineWidth(2);

  hgen_mupt->SetLineColor(4);
  hgen_mupt->SetLineWidth(2);

  //bo fsr
  hgen_elept_nofsr->SetLineColor(4);
  hgen_elept_nofsr->SetLineWidth(2);

  hgen_mupt_nofsr->SetLineColor(2);
  hgen_mupt_nofsr->SetLineWidth(2);

  hgenZee->SetLineColor(2);
  hgenZee->SetLineWidth(2);

  hgenZmm->SetLineColor(4);
  hgenZmm->SetLineWidth(2);

  hmva_zg_new->SetLineColor(4);
  hmva_zg_new->SetLineWidth(2);

  hmva_zg_old->SetLineColor(2);
  hmva_zg_old->SetLineWidth(2);
  
  //hmva_zg_new->Print();
  //hmva_zg_old->Print();
  hmva_zg_new->Scale(1./hmva_zg_new->Integral());
  hmva_zg_old->Scale(1./hmva_zg_old->Integral());

  TCanvas *c = new TCanvas("c", "c", 650, 500);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.01);
  pad1->Draw();
  pad1->cd();
  //pad1->SetLogy();

  hgen_elept->GetYaxis()->SetTitle("Entries/3 GeV");
  hgen_elept->GetXaxis()->SetTitle("p_{T}^{l} [GeV]");
  hgen_elept->Draw();
  hgen_mupt->Draw("same");
  TLegend *leg = new TLegend(0.6, 0.65, 0.8, 0.78);
  leg->AddEntry(hgen_elept, "electron", "f");
  leg->AddEntry(hgen_mupt, "muon", "f");
  leg->SetBorderSize(0);
  leg->Draw();

  TH1F *hratio = (TH1F*) hgen_elept->Clone();
  hratio->Divide(hgen_mupt);
  hratio->SetLineColor(1);

  c->cd();

  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1., 0.25);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();

  
  hratio->GetYaxis()->SetTitleOffset(1.6);
  hratio->SetMinimum(-6);
  hratio->SetTitleSize(0.14, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitle("ele/mu");
  hratio->GetXaxis()->SetTitle("p_{T}^{l} [GeV]");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0.8, 1.2);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->Draw();

  TLine *l = new TLine(-1,1,1,1);
  l->SetLineWidth(2);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");

  /*
  TCanvas *c2 = new TCanvas("c2", "c2", 650, 500);
  c2->cd();
  hgen_mupt->GetYaxis()->SetTitle("Entries/2 GeV");
  hgen_mupt->GetXaxis()->SetTitle("p_{T}^{l} [GeV]");
  hgen_mupt->Draw();
  hgen_mupt_nofsr->Draw("same");
  TLegend *leg1 = new TLegend(0.6, 0.6, 0.85, 0.72);
  leg1->AddEntry(hgen_mupt, "ZGToLLG_01J_5f", "f");
  leg1->AddEntry(hgen_mupt_nofsr, "ZGTo2LG", "f");
  leg1->SetBorderSize(0);
  leg1->Draw();
  */

  

  TCanvas *c1 = new TCanvas("c1", "c1", 650, 500);
  c1->cd();
  hmva_zg_new->GetYaxis()->SetTitle("Entries/0.1");
  hmva_zg_new->GetXaxis()->SetTitle("photon mva");
  hmva_zg_new->Draw();
  hmva_zg_old->Draw("same");
  //leg->Draw();


  TCanvas *c2 = new TCanvas("c2", "c2", 650, 500);
  c2->cd();
  hdral->GetYaxis()->SetTitle("Entries");
  hdral->GetXaxis()->SetTitle("#DeltaR(#gamma,l)");
  hdral->Draw();

  c3 = new TCanvas("c3", "c3", 650, 500);
  c3->cd();
  c3->SetLogy();
  hpho_fsr->GetYaxis()->SetTitle("Entries");
  hpho_fsr->GetXaxis()->SetTitle("p_{T}^{FSR} [GeV]");
  hpho_fsr->Draw();

  c4 = new TCanvas("c4", "c4", 650, 500);
  c4->cd();
  hdrlfsr->GetYaxis()->SetTitle("Entries");
  hdrlfsr->GetXaxis()->SetTitle("#DeltaR(FSR,lepton)");
  hdrlfsr->Draw();

  
  //c->SaveAs("plots/gen_subleadingPt_ele_mu_ZJets_AddAllFSR.pdf");
  c2->SaveAs("plots/deltaR_gamma_2ndEle_lhe_ZG.pdf");
  //c3->SaveAs("plots/FSR_pt_ZG.pdf");
  //c4->SaveAs("plots/dR_FSR_1stlepton_ZG.pdf");


}


