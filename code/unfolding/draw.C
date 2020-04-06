#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"


void draw(bool data = true, bool ele = true, bool barrel = true) {

  //from unfolding results
  TString infilename = "Zg_Unfolding_";
  if (data) infilename += "data_";
  else infilename += "MC_";
  if (ele) infilename += "ele_";
  else infilename += "mu_";
  if (barrel) infilename += "barrel_";
  else infilename += "endcap_";

  infilename += "Output_Under_OverFlow_SFs.root";

  TFile *inf = TFile::Open(infilename);

  TH1D *histUnfoldInput = (TH1D*) inf->Get("histUnfoldInput");
  TH1D *histUnfoldOutput = (TH1D*) inf->Get("histUnfoldOutput");
  TH1D *histFoldedBack = (TH1D*) inf->Get("histFolded");
  TH2D* hProbability = (TH2D*) inf->Get("hProbability");
  TH2D *histCorr = (TH2D*) inf->Get("histCorr");

  //from fitting
  int nbin = histUnfoldOutput->GetNbinsX();
  float yield_ele_EB[] = {9321.41,5163.81,2715.01,1521.63,1560.23,830.837,424.497,299.476,221.849,157.899,249.371,338.081};
  float yield_ele_EE[] = {3294.61,1566.42,800.275,500.93,571.496,267.454,163.537,106.199,72.4275,44.0274,84.5609,101.939};
  float yield_mu_EB[] = {20358.8,11726,6213.73,3071.73,3120.51,1627.88,970.789,603.005,440.653,288.724,414.577,600.015};
  float yield_mu_EE[] = {7713.54,3954.16,2059.68,1061.94,1144.7,555.72,367.064,213.475,150.939,108.463,134.131,160.611};

  float erryield_ele_EB[] = {134.684,89.9186,62.4727,45.7224,45.6447,31.8958,23.1947,19.228,15.5093,13.0264,16.3249,18.2567};
  float erryield_ele_EE[] = {103.507,62.3493,42.1064,29.1392,28.8518,19.5492,14.7705,11.8581,9.41864,7.03002,9.20017,10.0291};
  float erryield_mu_EB[] = {192.878,133.781,93.4323,66.9719,65.0807,45.3807,34.6132,26.9682,22.1344,17.7499,21.4794,24.4524};
  float erryield_mu_EE[] = {159.669,95.2252,62.2366,43.8492,41.5137,27.3889,21.4614,16.6414,12.9961,11.2009,12.3863,13.036};

  cout << "result of unfolding: " << endl;
  float ptbin[nbin+1];
  for (int i = 1; i <= nbin; i++) {
    ptbin[i-1] = histUnfoldOutput->GetXaxis()->GetBinLowEdge(i);
    cout << histUnfoldOutput->GetBinContent(i) << "\t";
  }
  ptbin[nbin] = histUnfoldOutput->GetXaxis()->GetBinUpEdge(nbin);

  cout << endl;

  cout << "error of unfolding: " << endl;
  for (int i = 1; i <= nbin; i++) {
    cout << histUnfoldOutput->GetBinError(i) << "\t";
  }

  cout << endl;

  TH1D *hinput_ele_EB = new TH1D("hinput_ele_EB", "hinput_ele_EB", nbin, ptbin);
  TH1D *hinput_ele_EE = new TH1D("hinput_ele_EE", "hinput_ele_EE", nbin, ptbin);
  TH1D *hinput_mu_EB = new TH1D("hinput_mu_EB", "hinput_mu_EB", nbin, ptbin);
  TH1D *hinput_mu_EE = new TH1D("hinput_mu_EE", "hinput_mu_EE", nbin, ptbin);

  for (int i = 1; i <= nbin; i++) {
    hinput_ele_EB->SetBinContent(i, yield_ele_EB[i-1]);
    hinput_ele_EE->SetBinContent(i, yield_ele_EE[i-1]);
    hinput_mu_EB->SetBinContent(i, yield_mu_EB[i-1]);
    hinput_mu_EE->SetBinContent(i, yield_mu_EE[i-1]);

    hinput_ele_EB->SetBinError(i, erryield_ele_EB[i-1]);
    hinput_ele_EE->SetBinError(i, erryield_ele_EE[i-1]);
    hinput_mu_EB->SetBinError(i, erryield_mu_EB[i-1]);
    hinput_mu_EE->SetBinError(i, erryield_mu_EE[i-1]);
  }

  //counting from MC
  TFile *finputMC = new TFile("ZgToLLg_5f_responseMatrix_MC_genbin.root", "read");
  TString hname;
  if (ele) {
    if (barrel) hname = "hrec_ele_EB";
    else hname ="hrec_ele_EE";
  }
  else {
    if (barrel) hname = "hrec_mu_EB";
    else hname ="hrec_mu_EE";
  }

  histUnfoldOutput->SetLineColor(2);
  histUnfoldOutput->SetLineWidth(2);

  TH1D *hinput;

  if (data) {
    if (ele && barrel) hinput = (TH1D*) hinput_ele_EB->Clone();
    else if (ele && !barrel) hinput = (TH1D*) hinput_ele_EE->Clone();
    else if (!ele && barrel) hinput = (TH1D*) hinput_mu_EB->Clone();
    else  hinput = (TH1D*) hinput_mu_EE->Clone();
  }
  else {
    hinput = (TH1D*) finputMC->Get(hname);
  }

  hinput->SetLineColor(4);
  hinput->SetLineWidth(2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetLogx();                                                                                                                                  
  pad1->SetFillColor(0);

  histUnfoldOutput->GetYaxis()->SetTitle("entries");
  histUnfoldOutput->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  histUnfoldOutput->GetXaxis()->SetLabelSize(0);
  histUnfoldOutput->Draw("e");
  //histFoldedBack->Draw("esames");
  hinput->Draw("esames"); 
  TLegend *lg = new TLegend(0.6, 0.7, 0.85, 0.82);
  lg->AddEntry(histUnfoldOutput, "unfolded", "fl");
  lg->AddEntry(hinput, "no unfolded", "fl");
  lg->SetBorderSize(0);
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->Draw();

  c1->cd();
  TH1F *hratio = (TH1F*) histUnfoldOutput->Clone();
  hratio->Divide(hinput);
  //hratio->SetLineColor(1);
  hratio->SetLineWidth(2);

  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  //pad2->SetLeftMargin(0.14);
  pad2->Draw();
  pad2->cd();
  pad2->SetLogx();
  pad2->SetGridy();

  hratio->GetYaxis()->SetTitleOffset(1.6);
  hratio->SetMinimum(-6);
  hratio->GetYaxis()->SetTitleOffset(1.1);
  hratio->SetTitleSize(0.1, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitle("red/blue");
  hratio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0.8, 1.2);
  hratio->GetYaxis()->SetNdivisions(505);
  hratio->GetXaxis()->SetTickLength(0.05);
  hratio->GetXaxis()->SetMoreLogLabels();

  hratio->Draw();

  setprecision(3);
  TCanvas *c2 = new TCanvas("c2", "c2", 650, 650);
  c2->cd();
  c2->SetLogy();
  c2->SetLogx();
  hProbability->GetXaxis()->SetTitleOffset(0.6);
  hProbability->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
  hProbability->GetYaxis()->SetTitle("p_{T}^{reco} [GeV]");
  hProbability->Rebin2D(1,2);
  hProbability->Draw("colz");


  TCanvas *c3 = new TCanvas("c3", "c3", 650, 650);
  c3->cd();
  c3->SetLogy();
  c3->SetLogx();
  //histCorr->GetZaxis()->SetRangeUser(0.1, 1);
  histCorr->Draw("colz");
  histCorr->Draw("box same");

  TString suf;
  if (data) suf += "_data_";
  else suf += "_MC_";
  if (ele) suf += "ele_";
  else suf += "mu_";
  if (barrel) suf += "barrel";
  else suf += "endcap";

  c1->SaveAs("plots/foldedresult" + suf + ".pdf");
  c2->SaveAs("plots/probabilities" + suf + ".pdf");
  c2->SaveAs("plots/correlation_coeff" + suf + ".pdf");



}
