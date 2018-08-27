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
//#include "CMS_lumi.C"

using namespace std;

void PhoMVAPlot(bool barrel = true) {

  gStyle->SetOptStat(0);

  //TFile *fmc = new TFile("../ana/minitrees_Zee/DYJets_amcatnlo_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr_BDTUpto6000.root", "read");
  //TFile *fmc = new TFile("../ana/minitrees_Zee/DYJets_amcatnlo_TnP_Zee_BDT_Upto6000_Corr_5Var_norho_nosieip.root", "read");
  TFile *fmc = new TFile("../ana/minitrees_Zee/DYJets_amcatnlo_TnP_Zee_BDT_Upto6000_5VarCorr.root", "read");
  TTree *tmc = (TTree*) fmc->Get("passingIdTree");

  TH1F *hphossmva = new TH1F("hphossmva", "", 20, -1, 1);
  TH1F *hphossmva_rw = new TH1F("hphossmva_rw", "", 20, -1, 1);
  TH1F *hdiffmva = new TH1F("hdiffmva", "", 50, -0.4, 0.8);
  TH2F *hmva_vs_rw = new TH2F("hmva_vs_rw", "mva vs mva ss correction", 20, -1, 1, 20, -1, 1);


  TCut cut = "passPhoID_Zg==1 && Zm>70 && Zm<110";
  if (barrel) cut += "fabs(Probe_SCEta)< 1.4442";
  else cut += "fabs(Probe_SCEta) > 1.566";

  tmc->Draw("Probe_ssmva:Probe_ssmva_rw >> hmva_vs_rw", cut);
  tmc->Draw("diff_mva >> hdiffmva", cut);

  cout << "plotting!!!!!!!!!!!" << endl;
  c1 = new TCanvas("c1", "c1", 650, 550);
  c1->cd();
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.15);
  hmva_vs_rw->GetYaxis()->SetTitle("MVA");
  hmva_vs_rw->GetXaxis()->SetTitle("showershape correction MVA");
  hmva_vs_rw->Draw("colz");

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(62);
  tx.SetTextSize(0.04);
  if (barrel) tx.DrawLatex(0.56, 0.88, "|#eta^{#gamma}| < 1.44");
  else tx.DrawLatex(0.56,0.88, "1.57 < |#eta^{#gamma}| < 2.5");

  c2 = new TCanvas("c2", "c2", 650, 550);
  c2->cd();
  hdiffmva->GetYaxis()->SetTitle("Events");
  hdiffmva->GetXaxis()->SetTitle("MVA - corrected MVA");
  hdiffmva->Draw("");
  if (barrel) tx.DrawLatex(0.62, 0.7, "|#eta^{#gamma}| < 1.4442");
  else tx.DrawLatex(0.62,0.7, "1.566 < |#eta^{#gamma}| < 2.5");
  tx.DrawLatex(0.62, 0.6, Form("mean = %.2f", hdiffmva->GetMean()));

  cout << "mean of histogram: " << hdiffmva->GetMean() << endl;

  //c1->SaveAs("plots/others/SSMVA_rw_vs_MVA_EB_5VarCorr.pdf");
  //c2->SaveAs("plots/others/diffMVA_EB_5VarCorr.pdf");


}


