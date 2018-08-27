#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCut.h>

#include <iostream>

using namespace std;

void histStyle(TH1F *h1) {
  h1->SetTitleSize(0.14, "XYZ");
  h1->SetLabelSize(0.14, "XYZ");
  h1->SetTitleFont(42, "XYZ");
  h1->GetYaxis()->SetTitle("data/ MC");
  //h1->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  h1->GetYaxis()->SetTitleOffset(0.55);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetRangeUser(0., 2);
  h1->GetYaxis()->SetNdivisions(5);
  //h1->Draw("ep");
}

void histStyle1 (TH1F *h) {
  h->SetFillColor(kCyan-7);
  h->SetLineColor(kCyan-7);
}

void histStyle2 (TH1F *h) {
  h->SetFillColor(kBlue);
  h->SetLineColor(kBlue);
}

void histStyle3 (TH1F *h) {
  h->SetFillColor(kBlue-9);
  h->SetLineColor(kBlue-9);
}

void histStyle4 (TH1F *h) {
  h->SetFillColor(kViolet+1);
  h->SetLineColor(kViolet+1);

}


void InputVarPlot(bool eb = true, TString hname = "hlepd0_lead") {

  gStyle->SetOptStat(0);

  cout << "get started" << endl;
  TString infname;
  if (eb) infname  = "histo/inputVar_Trainning/EB_";
  else infname  = "histo/inputVar_Trainning/EE_";

  const int nFile = 4;
  TFile *fin[nFile];
  TH1F *h[nFile];

  TH1F  *hMC = new TH1F();

  //fin[0] = new TFile(infname + "SingleEle_Run2016_TnP_Zee_chisocor_NewSSCor.root", "read");
  fin[0] = new TFile(infname + "SingleEle_Run2016_TnP_Zee_chisocor_NoGenMatch_NoSSCor.root", "read");
  //fin[1] = new TFile(infname + "DYJets_amcatnlo_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root", "read");
  //fin[2] = new TFile(infname + "TTbar_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root", "read");
  //fin[3] = new TFile(infname + "WJetsToLNu_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root", "read");

  fin[1] = new TFile(infname + "NoCor_DYJets_amcatnlo_TnP_Zee_BDT_Upto6000_Corr_EtaW_from_Sieie.root", "read");
  fin[2] = new TFile(infname + "NoCor_TTbar_TnP_Zee_BDT_Upto6000_Corr_EtaW_from_Sieie.root", "read");
  fin[3] = new TFile(infname + "NoCor_WJetsToLNu_TnP_Zee_BDT_Upto6000_Corr_EtaW_from_Sieie.root", "read");


  cout << "done opening files" << endl;
  for (int i = 0; i < nFile; i++) {
    h[i] = (TH1F*) fin[i]->Get(hname.Data());
  }

  for (int i = 1; i < nFile; i++) {
    if (i == 1) hMC = (TH1F*) h[i]->Clone();
    else hMC->Add(h[i]);
  }

  cout << "total entries from mc: " << hMC->Integral() << endl;
  cout << "total entries from data: " << h[0]->Integral() << endl;
  cout << "total entries from Zg: " << h[1]->Integral() << endl;
  cout << "total entries from Zj: " << h[2]->Integral() << endl;
  cout << "total entries from tt: " << h[3]->Integral() << endl;


  h[0]->Sumw2();
  TH1F *hRatio = (TH1F*) h[0]->Clone();
  hRatio->Divide(hMC);
  histStyle(hRatio);

  histStyle1(h[1]);
  histStyle2(h[2]);
  histStyle3(h[3]);

  THStack *hstack = new THStack("hstack", "");
  hstack->Add(h[3]);
  hstack->Add(h[2]);
  hstack->Add(h[1]);

  TString xname, yname, outname; 

  if (hname == "hphoHoE") {
    xname = "photon H/E";   yname = "Events"; outname = "plots/inputVarTrain/data_MC_Zee_phoHoE"; }
  if (hname == "hphosieie") {
    xname = "photon #sigma_{i#etai#eta}";   yname = "Events"; outname = "plots/inputVarTrain/data_MC_Zee_phoSieie"; }
  if (hname == "hphochiso") {
    xname = "charged iso^{#gamma}";   yname = "Events"; outname = "plots/inputVarTrain/data_MC_Zee_phoChargedIso"; }
  else if (hname == "hphowchiso") {
    xname = "worst charged iso^{#gamma}";   yname = "Events"; outname = "plots/inputVarTrain/data_MC_Zee_phoWorstChargedIso"; }
  if (hname == "hsieip") {
    xname = "photon #sigma_{i#etai#phi}"; yname = "Events"; outname = "plots/inputVarTrain/data_MC_Zee_phoSieip";}
  if (hname == "hetawidth") {
    xname = "#eta_{SC}^{width}"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_SCEtaWidth";}
  if (hname == "hphiwidth") {
    xname = "#phi_{SC}^{width}"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_SCPhiWidth";}
  if (hname == "hscRawE") {
    xname = "E^{raw}_{SC} [GeV]"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_SCRawE";}
  if (hname == "hs4") {
    xname = "E_{2x2}/E_{5x5}"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_S4";}
  if (hname == "hr9") {
    xname = "photon R9"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_phoR9";}
  if (hname == "hphophi") {
    xname = "#phi^{#gamma}"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_PhoPhi";}
  if (hname == "hphoeta") {
    xname = "#eta^{#gamma}"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_PhoSCEta";}
  if (hname == "hrho") {
    xname = "#rho"; yname = "Events"; outname ="plots/inputVarTrain/data_MC_Zee_rho";}
  if (hname == "hphossmva") {
    xname = "showershape MVA";   yname = "Events/0.1"; outname += "plots/inputVarTrain/data_MC_Zee_ssmva"; }
  

  //cout << "outname for saving: " << outname.Data() << endl;
  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.05);
  tx.SetTextFont(52);

  //plotting for canvas
  TCanvas* c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  c1->SetBottomMargin(0.15);
  c1->SetLogy();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  if (hname == "hphoHoE") pad1->SetLogy();
  h[0]->GetYaxis()->SetRangeUser(0.9, h[0]->GetMaximum()*2);
  h[0]->GetXaxis()->SetLabelSize(0);
  h[0]->GetYaxis()->SetTitle(yname);
  h[0]->Draw();
  hstack->Draw("histsame");
  h[0]->Draw("same");
  pad1->RedrawAxis();
  //TLegend *lg = new TLegend(0.44, 0.28, 0.7, 0.48);
  TLegend *lg = new TLegend(0.22, 0.56, 0.5, 0.7);
  lg->AddEntry(h[0], "data", "pe");
  lg->AddEntry(h[1], "DYJets", "f");
  lg->AddEntry(h[2], "T#barT", "f");
  lg->AddEntry(h[3], "WJets", "f");
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
  hRatio->GetXaxis()->SetTitle(xname);
  //hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hRatio->Draw("ep");
  tx.SetTextSize(0.08);
  tx.SetTextFont(42);
  tx.DrawLatex(0.5, 0.7, Form("#chi^{2}/NDF = %.2f", h[0]->Chi2Test(hMC, "UW CHI2/NDF")));


  TF1 *f1 = new TF1("myfuc", "pol0", -1, 1);
  //hRatio->Fit(f1);
  double ch2 = f1->GetChisquare();
  cout << "chi2: " << ch2 << endl;

  h[0]->Chi2Test(hMC, "UW P CHI2/NDF");
  cout << "chi2 between data and MC: " << h[0]->Chi2Test(hMC, "UW CHI2/NDF") << endl;

  if (eb)  c1->SaveAs(outname + "_EB_NoSSCorr.pdf");
  else c1->SaveAs(outname + "_EE_NoSSCorr.pdf");


}
