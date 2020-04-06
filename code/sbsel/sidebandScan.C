#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooAbsReal.h>
#include <RooSimultaneous.h>
#include <RooNDKeysPdf.h>
#include <RooKeysPdf.h>
#include <RooPlot.h>
#include "RooFitResult.h"

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace RooFit;
using namespace std;

void sidebandScan (bool ele = true, bool barrel = true, TString var = "pt", float minpt = 15., float maxpt = 20., float lower_sb = 3., float upper_sb = 5.) {

  TCut ptcut = Form("gamma_pt>%f && gamma_pt<%f", minpt, maxpt);
  TCut masscut = Form("boss_mass>%f && boss_mass<%f", minpt, maxpt);

  TCut cat;
  if (ele) cat = "leptType==11 && trig_Ele23_Ele12==1 && lept1_pt>20. && lept0_pt>25 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4";
  else cat = "leptType==13 && trig_Mu17_Mu8==1 && lept1_pt>20. && lept0_pt>25 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4";

  TCut cut;
  if (var == "pt") {
    if (barrel) cut = "isEB==1 && gamma_ChIso<2.";
    else cut ="isEE==1 && gamma_ChIso<1.5";

    cut += ptcut;
    //cut += cat;
  }
  
  else { //vs mass
    cut = "(isEB==1 && gamma_ChIso<2.) || (isEE==1 && gamma_ChIso<1.5)";
    cut += masscut; 
  }

  TCut cutbkg = Form("gamma_ChIso > %f && gamma_ChIso < %f", lower_sb, upper_sb);
  if (var == "pt") {
    if (barrel) cutbkg += "isEB==1";
    else cutbkg +="isEE==1";
    cutbkg += ptcut;
    //cutbkg += cat;
    cutbkg += "lept1_pt>20. && lept0_pt>25 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4 && ((leptType==11 && trig_Ele23_Ele12==1) || (leptType==13 && trig_Mu17_Mu8==1))";
  }
  else {
    cutbkg += masscut;
  }


  cout << "--- cut signal----:" ;
  cut.Print();
  cout << "--- cut background----:";
  cutbkg.Print();


  //TFile *fsig_templ = new TFile("../ana/minitrees/Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_cutatgen.root", "read");
  TFile *fsig_templ = new TFile("../ana/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  TFile *fbkg_templ = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_matching_dr0p2.root", "read");


  TTree *tsig_templ = (TTree*) fsig_templ->Get("outtree");
  TTree *tbkg_templ = (TTree*) fbkg_templ->Get("outtree");


  TH1F *hda = new TH1F("hda", "hist from data", 10, -1, 1);
  TH1F *hsig_templ = new TH1F("hsig_templ", "signal template", 10, -1, 1);
  TH1F *hbkg_templ = new TH1F("hbkg_templ", "background template", 10, -1, 1);
  TH1F *hbkg_templ_sb = new TH1F("hbkg_templ_sb", "background template in sb", 10, -1, 1);

  TH1F *hmc_bkg = new TH1F("hmc_bkg", "background", 10, -1, 1);
  TH1F *hmc_sig = new TH1F("hmc_sig", "signal", 10, -1, 1);

  tsig_templ->Draw("gamma_ssmva >> hsig_templ", cut+cat);
  tbkg_templ->Draw("gamma_ssmva >> hbkg_templ", cut+cat);
  tbkg_templ->Draw("gamma_ssmva >> hbkg_templ_sb", cutbkg);

  cout << "number of signal from template: " << hsig_templ->Integral() << endl;
  cout << "number of background from template side-band: " << hbkg_templ_sb->Integral() << endl;


  
  TString infda = "ExpectedYield/";
  if (var == "pt") infda += "Ptg/datacard/";
  else  infda += "Mllg/";
  infda += "data_";

  TString suf;
  if (var == "pt" ) {
    suf = Form("Pt%dto%d", (int)minpt, (int)maxpt);
    if (barrel) suf += "_EB";
    else suf += "_EE";
  }
  else if (var == "mass")
    suf = Form("Mllg%dto%d", (int)minpt, (int)maxpt);

  infda += suf;
  if (ele ) infda += "_ele.root";
  else infda += "_mu.root";

  //get expected signal and background 
  TFile *fin = TFile::Open(infda);
  TH1F *hsignal = (TH1F*) fin->Get("signal");
  TH1F *hda_obs = (TH1F*) fin->Get("data_obs");

  float inputSig = hsignal->Integral();
  float inputBkg = hda_obs->Integral() - hsignal->Integral();

  cout << "expected signal: " << inputSig << "t background: " << inputBkg << endl;


  /*********** close to signal (correct value) *************/
  /*
  if (var == "pt") {
    if (ele) {
      if ( barrel ) {
	if (minpt >= 15. && maxpt <= 20. ) { inputSig = 9992; inputBkg = 16921;}
	else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 5440; inputBkg = 6966;}
	else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 4339; inputBkg = 4821;}
	else if ( minpt >= 35. && maxpt <= 45. ) { inputSig = 1514; inputBkg = 1726;}
	else if ( minpt >= 45. && maxpt <= 65. ) { inputSig = 1331; inputBkg = 1129;}
	else if ( minpt >= 65. && maxpt <= 1000. ) { inputSig = 1254; inputBkg = 716;}
      }
      else {
	//if (minpt >= 15. && maxpt <= 20. ) { inputSig = 2841; inputBkg = 10248;}
	//if (minpt >= 15. && maxpt <= 20. ) { inputSig = 2144; inputBkg = 1838;} //test sieie < 0.027
	//if (minpt >= 15. && maxpt <= 20. ) { inputSig = 2399; inputBkg = 7084;} //charged iso < 0.5, H/E<0.023
	//if (minpt >= 15. && maxpt <= 20. ) { inputSig = 2132; inputBkg = 6205;} //charged iso < 0.034, H/E<0.023 (6205)
	if (minpt >= 15. && maxpt <= 20. ) { inputSig = 2551; inputBkg = 7122;} //sieie < 0.035
	//else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 1502; inputBkg = 3920;}
	else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 1448; inputBkg = 792;} //sieie < 0.027
	else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 1283; inputBkg = 2916;}
	else if ( minpt >= 35. && maxpt <= 45. ) { inputSig = 513; inputBkg = 993;}
	else if ( minpt >= 45. && maxpt <= 65. ) { inputSig = 439; inputBkg = 667;}
	else if ( minpt >= 65. && maxpt <= 1000. ) { inputSig = 489; inputBkg = 246;}

      }
    }

    else {// for muon channel
      if ( barrel ) {
	if (minpt >= 15. && maxpt <= 20. ) { inputSig = 23028; inputBkg = 37184;}
	else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 12497; inputBkg = 14730;}
	else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 9987; inputBkg = 9943;}
	else if ( minpt >= 35. && maxpt <= 45. ) { inputSig = 3329; inputBkg = 3405;}
	else if ( minpt >= 45. && maxpt <= 65. ) { inputSig = 2666; inputBkg = 2337;}
	else if ( minpt >= 65. && maxpt <= 1000. ) { inputSig = 2428; inputBkg = 1259;}
      }
      else {
	if (minpt >= 15. && maxpt <= 20. ) { inputSig = 6985; inputBkg = 22489;}
	else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 3897; inputBkg = 8569;}
	else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 3284; inputBkg = 5952;}
	else if ( minpt >= 35. && maxpt <= 45. ) { inputSig = 1106; inputBkg = 2137;}
	else if ( minpt >= 45. && maxpt <= 65. ) { inputSig = 918; inputBkg = 1359;}
	else if ( minpt >= 65. && maxpt <= 1000. ) { inputSig = 791; inputBkg = 609;}
      }
    }
  }

    //event as function of mass
  if (ele && var == "mass") {
    if (minpt >= 50. && maxpt <= 85. ) { inputSig = 1020; inputBkg = 149;}
    else if ( minpt >= 85. && maxpt <= 95. ) { inputSig = 12678; inputBkg = 900;} //bkg get from fit
    else if ( minpt >= 95. && maxpt <= 110. ) { inputSig = 2835; inputBkg = 5244;}
    else if ( minpt >= 110. && maxpt <= 135. ) { inputSig = 5602; inputBkg = 19765;}
    else if ( minpt >= 135. && maxpt <= 170. ) { inputSig = 4104; inputBkg = 12845;}
    else if ( minpt >= 170. && maxpt <= 210. ) { inputSig = 2120; inputBkg = 6435;}
    else if ( minpt >= 210. && maxpt <= 270. ) { inputSig = 1318; inputBkg = 4022;}
    else if ( minpt >= 270. && maxpt <= 350. ) { inputSig = 649; inputBkg = 1966;}
    else if ( minpt >= 350. && maxpt <= 1000. ) { inputSig = 600; inputBkg = 959;}
  }
  if (!ele && var == "mass") {
    if (minpt >= 50. && maxpt <= 85. ) { inputSig = 2093; inputBkg = 442;}
    else if ( minpt >= 85. && maxpt <= 95. ) { inputSig = 31846; inputBkg = 2000;} //bkg get from fit
    else if ( minpt >= 95. && maxpt <= 110. ) { inputSig = 6640; inputBkg = 10765;}
    else if ( minpt >= 110. && maxpt <= 135. ) { inputSig = 11999; inputBkg = 42564;}
    else if ( minpt >= 135. && maxpt <= 170. ) { inputSig = 8530; inputBkg = 27340;}
    else if ( minpt >= 170. && maxpt <= 210. ) { inputSig = 4073; inputBkg = 13494;}
    else if ( minpt >= 210. && maxpt <= 270. ) { inputSig = 2520; inputBkg = 8669;}
    else if ( minpt >= 270. && maxpt <= 350. ) { inputSig = 1239; inputBkg = 4090;}
    else if ( minpt >= 350. && maxpt <= 1000. ) { inputSig = 965; inputBkg = 1985;}
  }
  */

  
  
  //TH1F *hpull = new TH1F("hpull", "(fitted - generated)/generated (%)", 100, -20 ,20);
  TH1F *hpull = new TH1F("hpull", "pull of signal", 200, -20 ,20);
  //TH1F *hfitsig = new TH1F("hfitsig", "fitted signal", 500, 0.2*inputSig, 1.8*inputSig);
  TH1F *hfitsig = new TH1F("hfitsig", "fitted signal", 800, -2, 2);

  //generate toy MC
  int Ntoys = 1000;

  int count = 0;
  for (int i = 1 ; i<= Ntoys; i++) {
  TH1F *hsig_toy = new TH1F("hsig_toy", "", 10, -1, 1);
  hsig_toy->FillRandom(hsig_templ, inputSig);

  TH1F *hbkg_toy = new TH1F("hbkg_toy", "", 10, -1, 1);
  hbkg_toy->FillRandom(hbkg_templ, inputBkg);

  cout << "entries of toy signal: " << hsig_toy->GetEntries() << endl;
  cout << "entries of toy background: " << hbkg_toy->GetEntries() << endl;

  TH1F *hda_toy = (TH1F*) hsig_toy->Clone();
  hda_toy->Add(hbkg_toy);

  //cc = new TCanvas("cc", "cc", 650, 500);
  //cc->cd();
  //hsig_toy->Draw("e");

  //get data from histogram
  RooRealVar gamma_ssmva("gamma_ssmva", "photon shower shape mva", -1, 1);
  gamma_ssmva.setBins(10);

  //RooDataHist datah("datah", "# events in data", gamma_ssmva, Import(*hda));
  RooDataHist dahsig("dahsig", "# events in sigal", gamma_ssmva, Import(*hsig_templ));
  RooDataHist dahbkg("dahbkg", "# events in background", gamma_ssmva, Import(*hbkg_templ));

  RooDataHist datah("datah", "# events in data", gamma_ssmva, Import(*hda_toy));
  //RooDataHist dahsig("dahsig", "# events in sigal", gamma_ssmva, Import(*hsig_toy));
  //RooDataHist dahbkg("dahbkg", "# events in background", gamma_ssmva, Import(*hbkg_toy));

  //template pdf from histogram
  RooHistPdf sigpdf("sigpdf", "template pdf for signal", gamma_ssmva, dahsig);
  RooHistPdf bkgpdf("bkgpdf", "template pdf for background", gamma_ssmva, dahbkg);

  //background template from side-band
  RooDataHist dahbkg_sb("dahbkg_sb", "# events in background", gamma_ssmva, Import(*hbkg_templ_sb));
  RooHistPdf bkgpdf_sb("bkgpdf_sb", "template pdf for background in sb", gamma_ssmva, dahbkg_sb);


  RooRealVar frac("frac", "fraction of signal", 0.7);
  RooRealVar nsig("nsig", "no. signals", 200., 0., 200000.);
  RooRealVar nbkg("nbkg", "no. backgrounds", 200., 0., 200000);
  //RooRealVar nsig("nsig", "no. signals", inputSig);
  //RooRealVar nbkg("nbkg", "no. backgrounds", inputBkg);

  //model to data
  RooAddPdf model("model", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(nsig, nbkg));

  //model with bkg from sideband
  //RooRealVar nsig_sb("nsig_sb", "#signal event", 200., 0., 200000.);
  //RooRealVar nbkg_sb("nbkg_sb", "#background event", 200., 0., 200000.);

  float gensig = hsig_toy->GetEntries();
  float genbkg = hbkg_toy->GetEntries();

  RooRealVar nsig_sb("nsig_sb", "#signal event", gensig, 0.5*gensig, 1.5*gensig);
  RooRealVar nbkg_sb("nbkg_sb", "#background event", genbkg, 0.5*genbkg, 1.5*genbkg);
  //RooAddPdf model_sb("model_sb", "sig+bkg", RooArgList(sigpdf, bkgpdf_sb), RooArgList(nsig_sb, nbkg_sb));
  RooAddPdf model_sb("model_sb", "sig+bkg", RooArgList(bkgpdf_sb,sigpdf), RooArgList(nbkg_sb, nsig_sb));
  //RooAddPdf model_sb("model_sb", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(nsig_sb, nbkg_sb));

  //extended likelihood fit
  model.fitTo(datah, Extended());
  model_sb.fitTo(datah, Extended());

  RooPlot *frame = gamma_ssmva.frame(Title("fit data with signal and background template"));
  RooPlot *frame_sb = gamma_ssmva.frame(Title("fit data with signal and background template from sb"));


  cout << "generated signal: " << gensig << endl;
  cout << "fitted signal: " << nsig_sb.getVal() << " +/- " << nsig_sb.getError() << endl;
  cout << "fitted bkg: " << nbkg_sb.getVal() << " +/- " << nbkg_sb.getError() << endl;
  

  float pull_sig = (nsig_sb.getVal() - gensig) / nsig_sb.getError();
  //float pull_sig = (nsig.getVal() - gensig) / nsig.getError();
  //cout << "pull: " << pull_sig << endl;

  hpull->Fill(pull_sig);
  //hfitsig->Fill(nsig_sb.getVal());
  float bias = (nsig_sb.getVal() - gensig)/gensig;
  hfitsig->Fill(bias);


  if ( bias == 0.0) count++;

  datah.plotOn(frame, RooFit::Name("datah"), MarkerSize(1.0));
  model.plotOn(frame, RooFit::Name("model"), LineColor(kBlue), LineWidth(2));
  model.plotOn(frame, RooFit::Name("sigpdf"), Components("sigpdf"),LineColor(kRed), LineWidth(3), LineStyle(kDashed));
  model.plotOn(frame, RooFit::Name("bkgpdf"), Components("bkgpdf"),LineColor(kGreen), LineWidth(2));
  model.paramOn(frame, Format("NEU", AutoPrecision(2)), Layout(0.22, 0.77, 0.93));

  datah.plotOn(frame_sb, RooFit::Name("datah"), MarkerSize(1.0));
  model_sb.plotOn(frame_sb, RooFit::Name("model"), LineColor(kViolet), LineWidth(2));
  model_sb.plotOn(frame_sb, RooFit::Name("sigpdf"), Components("sigpdf"),LineColor(kRed), LineWidth(3), LineStyle(kDashed));
  model_sb.plotOn(frame_sb, RooFit::Name("bkgpdf"), Components("bkgpdf_sb"),LineColor(kGreen), LineWidth(2));
  model_sb.paramOn(frame_sb, Format("NEU", AutoPrecision(2)), Layout(0.22, 0.77, 0.93));


  if ( abs(bias) > 0.45 ) {
    cout << "bias: " << bias << endl;
  TCanvas *c1 = new TCanvas("c1", "c1", 1050, 550);
  c1->Divide(2,1);
  c1->cd(1);
  c1->SetLogy();
  frame->SetTitleOffset(1.5, "Y");
  frame->SetTitleOffset(1.3, "X");
  frame->GetYaxis()->SetRangeUser(0., frame->GetMaximum()*2);
  frame->Draw("same");
  TLegend *lg = new TLegend(0.3, 0.45, 0.56, 0.7);
  lg->AddEntry(frame->findObject("datah"),"toy MC","pe");
  lg->AddEntry(frame->findObject("model"),"model","l");
  lg->AddEntry(frame->findObject("sigpdf"),"signal","l");
  lg->AddEntry(frame->findObject("bkgpdf"),"background","l");
  lg->Draw();
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.04);
  tx.SetTextFont(42);
  tx.SetTextColor(kPink);
  if (barrel) tx.DrawLatex(0.3, 0.71, Form("EB, %d < p_{T}^{#gamma} < %d", (int)minpt, (int)maxpt));
  else tx.DrawLatex(0.3,0.71, Form("EE, %d < p_{T}^{#gamma} < %d", (int)minpt, (int)maxpt));
  tx.DrawLatex(0.3, 0.77, "bkg template from SR");

  c1->cd(2);
  frame_sb->SetTitleOffset(1.5, "Y");
  frame_sb->SetTitleOffset(1.3, "X");
  frame_sb->GetYaxis()->SetRangeUser(0., frame_sb->GetMaximum()*2);
  frame_sb->Draw("same");
  TLegend *lg_sb = new TLegend(0.3, 0.45, 0.56, 0.7);
  //lg_sb->AddEntry(frame_sb->findObject("toydata"),"pseudo-data","pe");
  lg_sb->AddEntry(frame_sb->findObject("datah"),"toyMC","pe");
  lg_sb->AddEntry(frame_sb->findObject("model"),"model","l");
  lg_sb->AddEntry(frame_sb->findObject("sigpdf"),"signal","l");
  lg_sb->AddEntry(frame_sb->findObject("bkgpdf"),"background","l");
  lg_sb->Draw();
  lg_sb->SetTextFont(42);
  lg_sb->SetTextSize(0.04);
  lg_sb->SetBorderSize(0);
  tx.DrawLatex(0.3, 0.74, "bkg template from SB");


  TString coutname = "fitting_toyMC";

  if (barrel) coutname += "_EB";
  else coutname += "_EE";

  coutname += Form("_PhoPt%dto%d", (int)minpt, (int)maxpt);
  coutname += Form("_chargedIso%dto%d", (int)lower_sb, (int)upper_sb);
  coutname += Form("_nBkg_%d_pull_%0.2f", (int) inputBkg, pull_sig);
  //if (ele) coutname += "_Ele";
  //else coutname += "_Mu";

  coutname += ".pdf";

  c1->SaveAs("plots/sibebandSel/fittingtoyMC/" + coutname);
  }

  }

  cout << "number of toy with bias < 5%: " << count << endl;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0111); 
  gStyle->SetStatY(0.70);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.18);
  gStyle->SetStatH(0.2);

  float lo_pull = hpull->GetXaxis()->GetBinLowEdge(hpull->FindFirstBinAbove(0.01,1));
  float hi_pull = hpull->GetXaxis()->GetBinUpEdge(hpull->FindLastBinAbove(0.01,1));


  cpull = new TCanvas("cpull", "cpull", 650, 650);
  cpull->cd();
  hpull->GetYaxis()->SetTitle("entries");
  hpull->GetXaxis()->SetTitle("pull of signal");
  hpull->GetXaxis()->SetRangeUser(lo_pull/2, hi_pull*1.5);
  hpull->Draw("e");

  TF1 *f1 = new TF1("gaus", "gaus", -4, 4);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  hpull->Fit(f1, "", "", lo_pull, hi_pull);

  //error on signal
  float lo_sig = hfitsig->GetXaxis()->GetBinLowEdge(hfitsig->FindFirstBinAbove(0.01,1));
  float hi_sig = hfitsig->GetXaxis()->GetBinUpEdge(hfitsig->FindLastBinAbove(0.01,1));


  cfitsig = new TCanvas("cfitsig", "cfitsig", 650, 650);
  cfitsig->cd();
  hfitsig->GetYaxis()->SetTitle("entries");
  hfitsig->GetXaxis()->SetTitle(Form("fitted signal (input %d)", (int) inputSig));
  hfitsig->GetXaxis()->SetRangeUser(lo_sig/1.5, hi_sig*1.2);
  //hfitsig->GetXaxis()->SetTitle("#frac{fit-gen}{gen}");
  //hfitsig->GetXaxis()->SetRangeUser(-0.5, 0.5);
  hfitsig->Draw("e");

  TF1 *f2 = new TF1("gaus", "gaus", lo_sig, hi_sig);
  f2->SetLineColor(kBlue);
  f2->SetLineWidth(2);
  //hfitsig->Fit(f2, "", "", lo_sig, hi_sig);
  hfitsig->Fit(f2, "", "", lo_sig, hi_sig);

  cout << "================= result ===================" << endl;
  cout << "mean of pull signal: " << setprecision(2) << f1->GetParameter(1) << endl;
  cout << "fitted signal: " << setprecision(4) << f2->GetParameter(1) << endl;
  cout << "input of signal: " << setprecision(4) << inputSig << endl;
  //cout << "error of template: " << setprecision(3) << abs(f2->GetParameter(1) - inputSig)/inputSig *100 << endl;
  cout << "error of template: " << setprecision(3) << f2->GetParameter(1)*100 << endl;

  ofstream fout;
  string filename;
  filename = "biasValues/pull_error_";
  filename += Form("PhoPt_%dto%d", (int) minpt, (int) maxpt);
  if (barrel) filename += "_EB";
  else filename += "_EE";
  //filename += "_QCD50to80_Signal_Zg";
  filename += "_ZJets_Signal_Zg";
  if (ele) filename += "_BDT_up6000_5var_ele_newZG.txt";
  else filename += "_BDT_up6000_5varmu_newZG.txt";
  fout.open(filename, ios::app);

  //fout << (int)lower_sb << " < iso < " << (int)upper_sb << "\t " << setprecision(2) << f1->GetParameter(1) 
  // << "\t \t " << abs(f2->GetParameter(1) - inputSig)/inputSig *100 << endl;
  fout << (int)lower_sb << " < iso < " << (int)upper_sb << "\t " << setprecision(2) << f1->GetParameter(1) 
       << "\t \t " << f2->GetParameter(1)*100 << endl;


  //for saving
  //TString dir = "plots/sibebandSel/SignalPull/selfbias/";
  TString dir = "plots/sibebandSel/BDT_up6000/bias/";
  //TString subdir = Form("PhoPt%dto%d/QCD50to80/muon/", (int)minpt, (int)maxpt);
  TString subdir = Form("PhoPt%dto%d/ZJets_newZG/ele/", (int)minpt, (int)maxpt);
  //TString subdir = "Unc/";
  const int dir_err = system("mkdir -p " + dir + subdir);
  if (-1 == dir_err)
    {
      printf("Error creating directory!n");
      exit(1);
    }



  TString outname = "";
  if (barrel) outname += "_EB";
  else outname += "_EE";

  outname += Form("_PhoPt%dto%d", (int)minpt, (int)maxpt);
  outname += Form("_chargedIso%dto%d", (int)lower_sb, (int)upper_sb);
  //outname += Form("_nBkg_%d", (int) inputBkg);
  if (ele) outname += "_Ele";
  else outname += "_Mu";

  outname += ".pdf";


  cpull->SaveAs(dir + subdir + "bias" + outname);
  cfitsig->SaveAs(dir + subdir + "FittedSignal" + outname);
  


}
