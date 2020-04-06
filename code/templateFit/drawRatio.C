#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

using namespace std;

void drawRatio(bool ele = true, bool barrel = true, TString cat= "pt", float minpt = 15, float maxpt = 20, bool FirstFit = true) {

  gStyle->SetOptStat(0);

  TString suf;
  if (cat == "pt" ) {
    //suf = Form("Pt%dto%d", (int)minpt, (int)maxpt);
    suf = Form("Ptllg%dto%d", (int)minpt, (int)maxpt);
  } 
  else if (cat == "mass") 
    suf = Form("Mllg%dto%d", (int)minpt, (int)maxpt);

  if (barrel) suf += "_EB";
  else suf += "_EE";

  TString infda, in_dir;

  //infda = "SB_EB7to13_EE6to14/";
  //in_dir = "SB_EB7to13_EE6to14/";
  infda = "SB_EB7to13_EE6to14/exclusive/";
  in_dir = "SB_EB7to13_EE6to14/exclusive/";
  //infda = "SB_EB7to13_EE5to11/";
  //in_dir = "SB_EB7to13_EE5to11/";

  //if (cat == "pt") infda += "Ptg/datacard/";
  if (cat == "pt") infda += "Ptllg/datacard/";
  else  infda += "Mllg/datacard/";
  infda += "data_";

  infda += suf;
  if (ele ) infda += "_ele.root";
  else infda += "_mu.root";


  //fitting file
  TString infit ;
  if (FirstFit ) {
    //if ( cat== "pt") infit = "SB_EB7to13_EE6to14/Ptg/1stFit/fitDiagnostics_";
    if ( cat== "pt") infit = "SB_EB7to13_EE6to14/Ptllg/1stFit/fitDiagnostics_";
    else infit = "SB_EB7to13_EE6to14/Mllg/1stFit/fitDiagnostics_";
  }

  else { 
    //if ( cat== "pt") infit = in_dir + "/Ptg/2ndFit/fitting/fitDiagnostics_";
    if ( cat== "pt") infit = in_dir + "/Ptllg/2ndFit/fitting/fitDiagnostics_";
    else infit = in_dir + "/Mllg/2ndFit/fitting/fitDiagnostics_";
  }

  infit += suf;
  if (ele ) infit += "_ele.root";
  else infit += "_mu.root";

  cout << "input fitting file: " << infit << endl;
  TFile *fda = new TFile(infda.Data(), "read");
  TFile *temp_fit = new TFile(infit.Data(), "read"); 

  cout << "input file of data: " << infda.Data() << endl;
  cout << "input fitting file: " << infit.Data() << endl;

  TH1F *hda = (TH1F*) fda->Get("data_obs");

  TH1D *hsignal_fit = new TH1D ("hsignal_fit", "fitted signal", 10, -1, 1);
  TH1D *hbackground_fit = new TH1D ("hbackground_fit", "fitted background", 10, -1, 1);
  TH1D *htotal_fit = new TH1D ("htotal_fit", "fitted data", 10, -1, 1);

  TH1D *hsignal_fit_temp = (TH1D*) temp_fit->Get("shapes_fit_s/total_signal");
  TH1D *hbackground_fit_temp = (TH1D*) temp_fit->Get("shapes_fit_s/total_background");
  TH1D *htotal_fit_temp = (TH1D*) temp_fit->Get("shapes_fit_s/total_overall");

  //cout << "get signal yield after first fit" << endl;
  //event yield from fit_Result
  RooArgSet *norm_fit_s = (RooArgSet*) temp_fit->Get("norm_fit_s");
  RooRealVar* parfit = (RooRealVar*) norm_fit_s->find("Zg/signal");
  float  sigYield = parfit->getValV();
  float sigYield_err = parfit->getError();

  float sig_event = hsignal_fit_temp->Integral();
  float rescale = 1.;
  if (sigYield > 0.) rescale = sigYield/sig_event;

  float ymax = 0;

  //hsignal_fit_temp->Print("all");

  for (int i=1; i<=htotal_fit->GetNbinsX(); i++) {
    if((hda->GetBinContent(i)+hda->GetBinError(i))>ymax) ymax = hda->GetBinContent(i)+hda->GetBinError(i);
    hsignal_fit->SetBinContent(i, hsignal_fit_temp->GetBinContent(i)*rescale);
    hbackground_fit->SetBinContent(i, hbackground_fit_temp->GetBinContent(i)*rescale);
    htotal_fit->SetBinContent(i, htotal_fit_temp->GetBinContent(i)*rescale);

    //cout << "error of fit: " << hsignal_fit->GetBinError(i) << endl;
    hsignal_fit->SetBinError(i, hsignal_fit_temp->GetBinError(i)*rescale);
    hbackground_fit->SetBinError(i, hbackground_fit_temp->GetBinError(i)*rescale);
    htotal_fit->SetBinError(i, htotal_fit_temp->GetBinError(i)*rescale);

  }

  //hbackground_fit->Print("all");

  Double_t chi2=0.;
  //chi2=hda->Chi2Test(htotal_fit,"P CHI2");
  chi2=hda->Chi2Test(htotal_fit,"CHI2/NDF");


  hda->Sumw2();
  hda->SetMarkerSize(1);
  hda->SetMarkerStyle(20);

  htotal_fit->SetLineColor(kBlue);
  htotal_fit->SetLineWidth(2);

  hsignal_fit->SetLineColor(kCyan-3);
  hsignal_fit->SetFillColor(kCyan-3);
  hsignal_fit->SetLineWidth(1);

  hbackground_fit->SetLineColor(kRed-4);
  hbackground_fit->SetLineWidth(1);
  hbackground_fit->SetFillColor(kRed-4);
  hbackground_fit->SetFillStyle(3002);


  gStyle->SetOptStat(0);

  cout << "plotting canvas" << endl;
  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  //pad1->SetLogy();
  pad1->SetFillColor(0);

  htotal_fit->GetYaxis()->SetTitleOffset(1.6);
  htotal_fit->SetTitleSize(0.05, "XYZ");
  htotal_fit->SetLabelSize(0.05, "XYZ");
  htotal_fit->SetTitleFont(42, "XYZ");
  htotal_fit->GetYaxis()->SetTitle("Entries/0.1");

  htotal_fit->SetLabelSize(0, "X");
  htotal_fit->SetTitleOffset(1.5, "Y");
  htotal_fit->SetTitleOffset(1.3, "X");
  htotal_fit->GetYaxis()->SetRangeUser(1, htotal_fit->GetMaximum()*15);
  htotal_fit->SetMaximum(ymax*1.2);
  //if (maxpt >= 75. && barrel) htotal_fit->SetMaximum(ymax*1.8);
  if (maxpt >= 75. && barrel) htotal_fit->GetYaxis()->SetRangeUser(0, htotal_fit->GetMaximum()*1.8);
  if (cat == "pt" && minpt >= 120) htotal_fit->GetYaxis()->SetRangeUser(0, htotal_fit->GetMaximum()*1.8);
  htotal_fit->Draw("hist");
  hsignal_fit->Draw("hist same");
  hbackground_fit->Draw("hist same");
  htotal_fit->Draw("hist same");
  hda->Draw("same");
  TLegend *lg;
  if (cat == "pt") lg = new TLegend(0.4, 0.4, 0.8, 0.64);
  if (cat == "pt" && minpt >= 120) lg = new TLegend(0.28, 0.4, 0.8, 0.64);
  else lg = new TLegend(0.28, 0.5, 0.68, 0.74);
  lg->AddEntry(hda,Form("data    %.1f events", hda->Integral()),"pe");
  lg->AddEntry(htotal_fit,Form("fitted   %.1f events", htotal_fit->Integral()),"l");
  //lg->AddEntry(hsignal_fit,Form("signal  %.1f #pm %.1f events", hsignal_fit->Integral(), sig_err),"f");
  lg->AddEntry(hsignal_fit,Form("signal  %.1f #pm %.1f events", sigYield, sigYield_err),"f");
  lg->AddEntry(hbackground_fit,"Background","f");
  lg->Draw();
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  
  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.04);
  tx.SetTextFont(42);
  if (cat == "pt" ) tx.DrawLatex(0.4, 0.35, Form("#chi^{2}/NDF = %.2f", chi2));  
  else tx.DrawLatex(0.34, 0.45, Form("#chi^{2}/NDF = %.2f", chi2));

  tx.SetTextSize(0.042);
  if (ele) tx.DrawLatex(0.42, 0.88, "Z#gamma#rightarrow ee#gamma");
  else tx.DrawLatex(0.4, 0.88, "Z#gamma #rightarrow #mu#mu#gamma");
  //if (cat == "pt" ) {
    if (barrel) tx.DrawLatex(0.4, 0.82, "0 < |#eta^{#gamma}| < 1.5");
    else tx.DrawLatex(0.4,0.82, "1.5 < |#eta^{#gamma}| < 2.5");
    //if (cat == "pt" ) tx.DrawLatex(0.4,0.76, Form("%d < p_{T}^{#gamma} < %d", (int)minpt, (int)maxpt));
    if (cat == "pt" ) tx.DrawLatex(0.4,0.76, Form("%d < p_{T}^{ll#gamma} < %d", (int)minpt, (int)maxpt));
    //}
  else tx.DrawLatex(0.4,0.76, Form("%d < M_{ll#gamma} < %d", (int)minpt, (int)maxpt));

  //cout << "ratio between fitting and data" << endl;
  c1->cd();
  TH1F *hratio = (TH1F*) hda->Clone();
  hratio->Divide(htotal_fit);
  hratio->SetLineColor(1);

  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();

  hratio->GetYaxis()->SetTitleOffset(1.6);
  hratio->SetMinimum(-6);
  hratio->SetTitleSize(0.14, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitle("data/ fit");
  hratio->GetXaxis()->SetTitle("photon MVA");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0.5, 1.5);
  hratio->GetYaxis()->SetNdivisions(8);
  hratio->Draw();

  TLine *l = new TLine(-1,1,1,1);
  l->SetLineWidth(2);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");


  //----------------------------------------------------------//
  //--------------- bias correction -------------------------//


  string chan;
  
  //for inclusive
  if (cat == "pt") {
    //if (ele && barrel) chan = "mean_width_PhoPt20_Ptg_EB_SB_7to13_ele.txt";
    //else if (ele && !barrel) chan = "mean_width_PhoPt20_Ptg_EE_SB_6to14_ele.txt";
    //else if (!ele && barrel) chan = "mean_width_PhoPt20_Ptg_EB_SB_7to13_mu.txt";
    //else if (!ele && !barrel) chan = "mean_width_PhoPt20_Ptg_EE_SB_6to14_mu.txt";

    if (ele && barrel) chan = "mean_width_PhoPt20_Ptllg_EB_SB_7to13_excl_ele.txt";
    else if (ele && !barrel) chan = "mean_width_PhoPt20_Ptllg_EE_SB_6to14_excl_ele.txt";
    else if (!ele && barrel) chan = "mean_width_PhoPt20_Ptllg_EB_SB_7to13_excl_mu.txt";
    else if (!ele && !barrel) chan = "mean_width_PhoPt20_Ptllg_EE_SB_6to14_excl_mu.txt";
  }
  else if (cat == "mass") {
    if (ele && barrel) chan = "mean_width_PhoPt20_Mllg_EB_SB_7to13_ele_rebin.txt";
    else if (ele && !barrel) chan = "mean_width_PhoPt20_Mllg_EE_SB_6to14_ele_rebin.txt";
    else if (!ele && barrel) chan = "mean_width_PhoPt20_Mllg_EB_SB_7to13_mu_rebin.txt";
    else if (!ele && !barrel) chan = "mean_width_PhoPt20_Mllg_EE_SB_6to14_mu_rebin.txt";
  }


  /*
  //for exclusive
  if (cat == "pt") {
    if (ele && barrel) chan = "mean_width_PhoPt20_Ptg_EB_SB_7to13_excl_ele.txt";
    else if (ele && !barrel) chan = "mean_width_PhoPt20_Ptg_EE_SB_6to14_excl_ele.txt";
    else if (!ele && barrel) chan = "mean_width_PhoPt20_Ptg_EB_SB_7to13_excl_mu.txt";
    else if (!ele && !barrel) chan = "mean_width_PhoPt20_Ptg_EE_SB_6to14_excl_mu.txt";
  }
  else if (cat == "mass") {
    if (ele && barrel) chan = "mean_width_PhoPt20_Mllg_EB_SB_7to13_excl_ele.txt";
    else if (ele && !barrel) chan = "mean_width_PhoPt20_Mllg_EE_SB_6to14_excl_ele.txt";
    else if (!ele && barrel) chan = "mean_width_PhoPt20_Mllg_EB_SB_7to13_excl_mu.txt";
    else if (!ele && !barrel) chan = "mean_width_PhoPt20_Mllg_EE_SB_6to14_excl_mu.txt";
  }
  */

  ifstream ferr;
  ferr.open(chan.c_str());
  cout << "file of mean and width  " << chan << endl;

  if (!ferr) {
    cout << "file could not be opened" << endl;
    exit;
  }

  float mean = 0., width = 0.;
  string out;
  while (getline(ferr, out)) {
    if (out.find("mean") != string::npos) {
      stringstream ss(out);
      string s1, s2, s3, s4, s5;
      int a1, a2;
      float _mean, _w;
      ss >> a1 >> s1 >> s2  >> s3 >> a2 >> s4 >> _mean >> s5 >> _w;
      if ((int)minpt == a1 && (int)maxpt == a2) {mean = _mean; width = _w;}
      cout << _mean << endl;
    }
  }

  cout << "bias corrected signal yield: " << sigYield/(1+mean) << endl;
  //---------------------done bias corrected yield-----------------------//


  ofstream fout;
  string filename;
  filename = "results/signal_yield";
  //if (cat == "pt" ) filename += "_Ptg";
  if (cat == "pt" ) filename += "_Ptllg";
  else filename += "_Mllg";

  if (barrel) filename += "_EB_SB_7to13";
  else filename += "_EE_SB_6to14";

  //if (ele) filename += "_newZG_newSB_ele_TUnfold_bin_isPVGood.txt";
  //else filename += "_newZG_newSB_mu_TUnfold_bin_isPVGood.txt";

  if (ele) filename += "_ele_excl_biasDn.txt";
  else filename += "_mu_excl_biasDn.txt";

  fout.open(filename, ios::app);

  mean = mean-width;
  //cout << "mean: " << mean  << endl;
  if (cat == "pt" ) {
    float tot = htotal_fit->Integral();
    fout << (int) minpt << " < pt < " << (int) maxpt << "\t " << sigYield << " +/- " << sigYield_err 
    	 << "\t bias corrected: " << sigYield/(1+mean) << " +/- " << sigYield_err/(1+mean) << "\t fakerate: " << sigYield/((1+mean)*tot) 
	 << " +/- " << sigYield_err/((1+mean)*tot) << endl;

    //float tot = htotal_fit->Integral();
    //fout << (int) minpt << " < pt < " << (int) maxpt << "\t " << sigYield << " +/- " << sigYield_err << "\t fakerate: " << sigYield/tot 
    // << " +/- " << sigYield_err/tot << endl;
  }
  else
    fout << (int) minpt << " < mllg < " << (int) maxpt << "\t " << sigYield << " +/- " << sigYield_err 
    	 << "\t bias corrected: " << sigYield/(1+mean) << " +/- " << sigYield_err/(1+mean) << endl;

  //cout << "saving plot" << endl;

  //TString dir = "plots/fitting/EB7to13_EE6to14/Refit_April11/";
  //TString dir = "plots/fitting/EB7to13_EE6to14/Refit_April24/"; //for inclusive
  //TString dir = "plots/fitting/EB7to13_EE6to14/exclusive/"; //for inclusive
  TString dir = "plots/fitting/EB7to13_EE6to14/";

  //if ( cat == "pt") dir += "Ptg/";
  if ( cat == "pt") dir += "Ptllg/";
  else dir+= "Mllg/";

  const int dir_err = system("mkdir -p " + dir);
  if (-1 == dir_err)
    {
      printf("Error creating directory!n");
      exit(1);
    }

  TString outname = "fitting_data";
  if (cat == "pt" ) {
    //outname += Form("_PhoPt%dto%d", (int)minpt, (int)maxpt);
    outname += Form("_PhoPtllh%dto%d", (int)minpt, (int)maxpt);
    //if (barrel) outname += "_EB";
    //else outname += "_EE";
  }
  else 
    outname += Form("_Mllg%dto%d", (int)minpt, (int)maxpt);
  if (barrel) outname += "_EB";
  else outname += "_EE";

  if (ele) outname += "_ele";
  else outname += "_mu";


  if (FirstFit)  outname += "_newZG_NoLepFlavorForSB_1stFit.pdf";
  else outname += "_2ndFit_newZG_NoLepFlavorForSB_excl.pdf";

  //c1->SaveAs(dir + outname);



}
