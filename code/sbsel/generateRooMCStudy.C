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

using namespace RooFit;
using namespace std;

void generateRooMCStudy (bool ele = true, bool barrel = true, float minpt = 15., float maxpt = 20., float lower_sb = 3., float upper_sb = 5.) {

  TCut ptcut = Form("gamma_pt>%f && gamma_pt<%f", minpt, maxpt);

  TCut cut;
  if (barrel) cut = "isEB==1 && gamma_ChIso<2.";
  else cut ="isEE==1 && gamma_ChIso<1.5";
  cut += ptcut;

  TCut cutbkg = Form("gamma_ChIso > %f && gamma_ChIso < %f", lower_sb, upper_sb);
  if (barrel) cutbkg += "isEB==1";
  else cutbkg +="isEE==1";
  cutbkg += ptcut;

  TCut cat;
  if (ele) cat = "leptType==11";
  else cat = "leptType==13";

  //cut += cat;

  cout << "--- cut signal----:" ;
  cut.Print();
  cout << "--- cut background----:";
  cutbkg.Print();


  TFile *fsig_templ = new TFile("output/bkgTemplate_summer16_420_job_summer16_WGToLNuG_amcatnlo_ext3.root", "read");
  TFile *fbkg_templ = new TFile("output/bkgTemplate_summer16_420_job_summer16_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root", "read");
  //TFile *fbkg_templ = new TFile("output/bkgTemplate_summer16_420_job_summer16_QCD_Pt_80to120_TuneCUETP8M1_pythia8.root", "read");
  //TFile *fbkg_templ = new TFile("output/bkgTemplate_summer16_420_LooseSieie_job_summer16_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root", "read");
  //TFile *fbkg_templ = new TFile("output/bkgTemplate_summer16_420_LooseSieie_job_summer16_QCD_Pt_80to120_TuneCUETP8M1_pythia8.root", "read");


  TFile *fmc_sig = new TFile("output/Zg_aMCatNLO_SortedPt_geniso.root", "read");
  TFile *fmc_bkg = new TFile("output/ZJets_aMCatNLO_SortedPt.root", "read");
  TChain *tda = new TChain("data");

  tda->Add("output/Zg_aMCatNLO_SortedPt_geniso.root/outtree");
  tda->Add("output/ZJets_aMCatNLO_SortedPt.root/outtree");


  TTree *tsig_templ = (TTree*) fsig_templ->Get("outtree");
  TTree *tbkg_templ = (TTree*) fbkg_templ->Get("outtree");

  //for checking
  TTree *tmc_bkg = (TTree*) fmc_bkg->Get("outtree");
  TTree *tmc_sig = (TTree*) fmc_sig->Get("outtree");

  TH1F *hda = new TH1F("hda", "hist from data", 20, -1, 1);
  TH1F *hsig_templ = new TH1F("hsig_templ", "signal template", 20, -1, 1);
  TH1F *hbkg_templ = new TH1F("hbkg_templ", "background template", 20, -1, 1);
  TH1F *hbkg_templ_sb = new TH1F("hbkg_templ_sb", "background template in sb", 20, -1, 1);

  TH1F *hmc_bkg = new TH1F("hmc_bkg", "background", 20, -1, 1);
  TH1F *hmc_sig = new TH1F("hmc_sig", "signal", 20, -1, 1);

  //tda->Draw("gamma_ssmva >> hda", cut+cat);
  tda->Draw("gamma_ssmva >> hda", cut );
  tsig_templ->Draw("gamma_ssmva >> hsig_templ", cut);
  tbkg_templ->Draw("gamma_ssmva >> hbkg_templ", cut);
  tbkg_templ->Draw("gamma_ssmva >> hbkg_templ_sb", cutbkg);

  //check distribution
  tmc_bkg->Draw("gamma_ssmva >> hmc_bkg", cut);
  tmc_sig->Draw("gamma_ssmva >> hmc_sig", cut);

  float inputSig = 200., inputBkg = 200.;

  if ( barrel ) {
    if (minpt >= 15. && maxpt <= 20. ) { inputSig = 7156.; inputBkg = 13219.;}
    //if (minpt >= 15. && maxpt <= 20. ) { inputSig = 7156.; inputBkg = 3000.;}
    else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 4032.; inputBkg = 5426.;}
    else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 3434.; inputBkg = 3690.;}
    else if ( minpt >= 35. && maxpt <= 45. ) { inputSig = 1401.; inputBkg = 1173.;}
    else if ( minpt >= 45. && maxpt <= 65. ) { inputSig = 1350.; inputBkg = 604.;}
    else if ( minpt >= 65. && maxpt <= 1000. ) { inputSig = 1361.; inputBkg = 294.;}
  }
  else {
    if (minpt >= 15. && maxpt <= 20. ) { inputSig = 1563.; inputBkg = 2345.;}
    else if ( minpt >= 20. && maxpt <= 25. ) { inputSig = 1112.; inputBkg = 1050.;}
    else if ( minpt >= 25. && maxpt <= 35. ) { inputSig = 1054.; inputBkg = 897.;}
    else if ( minpt >= 35. && maxpt <= 55. ) { inputSig = 873.; inputBkg = 366.;}
    else if ( minpt >= 55. && maxpt <= 1000. ) { inputSig = 812.; inputBkg = 366.;}
    cout << "----input of signal in endcaps-----: " << inputSig << endl;
  }

  //inputSig = 26000., inputBkg = 14000.;
  cout << "input of signal: " << inputSig << endl;

  //get data from histogram
  RooRealVar gamma_ssmva("gamma_ssmva", "photon shower shape mva", -1, 1);
  gamma_ssmva.setBins(20);

  RooDataHist datah("datah", "# events in data", gamma_ssmva, Import(*hda));
  RooDataHist dahsig("dahsig", "# events in sigal", gamma_ssmva, Import(*hsig_templ));
  RooDataHist dahbkg("dahbkg", "# events in background", gamma_ssmva, Import(*hbkg_templ));

  //template pdf from histogram
  RooHistPdf sigpdf("sigpdf", "template pdf for signal", gamma_ssmva, dahsig);
  RooHistPdf bkgpdf("bkgpdf", "template pdf for background", gamma_ssmva, dahbkg);

  //background template from side-band
  RooDataHist dahbkg_sb("dahbkg_sb", "# events in background", gamma_ssmva, Import(*hbkg_templ_sb));
  RooHistPdf bkgpdf_sb("bkgpdf_sb", "template pdf for background in sb", gamma_ssmva, dahbkg_sb);


  //double frac = 0.4;

  RooRealVar frac("frac", "fraction of signal", 0.7);
  //RooRealVar nsig("nsig", "no. signals", 200., 0., 200000.);
  //RooRealVar nbkg("nbkg", "no. backgrounds", 400., 0., 200000);
  RooRealVar nsig("nsig", "no. signals", inputSig);
  RooRealVar nbkg("nbkg", "no. backgrounds", inputBkg);

  //model to data
  RooAddPdf model("model", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(nsig, nbkg));
  //RooAddPdf model("model", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(frac));

  RooAddPdf model_sig("model_sig", "sig", RooArgList(sigpdf), RooArgList(nsig));
  RooAddPdf model_bkg("model_bkg", "bkg", RooArgList(bkgpdf), RooArgList(nbkg));

  //model with bkg from sideband
  RooRealVar nsig_sb("nsig_sb", "#signal event", 400., 0., 200000.);
  RooRealVar nbkg_sb("nbkg_sb", "#background event", 200., 0., 200000.);
  //RooAddPdf model_sb("model_sb", "sig+bkg", RooArgList(sigpdf, bkgpdf_sb), RooArgList(nsig_sb, nbkg_sb));
  RooAddPdf model_sb("model_sb", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(nsig_sb, nbkg_sb));

  //generate toy data
  //RooDataHist *toydata = model.generateBinned(gamma_ssmva, 40000);

  
  //RooMCStudy *toydata = new RooMCStudy(model,gamma_ssmva,Binned(kTRUE),Silence(),Extended(), FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  //RooMCStudy *toydata = new RooMCStudy(model,gamma_ssmva,FitModel(model_sb),Silence(),Binned(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  RooMCStudy *toydata = new RooMCStudy(model,gamma_ssmva,FitModel(model_sb),Silence(),Binned(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
  RooMCStudy *toydata1 = new RooMCStudy(model,gamma_ssmva,FitModel(model_sb),Silence(),Binned(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));


  UInt_t Ntoys = 1000;
  //toydata->generateAndFit(Ntoys, 1000, kTRUE, "generatedata/toymc_%04d.dat");
  toydata->generateAndFit(Ntoys, int(inputSig+inputBkg), kTRUE);
  //toydata1->generateAndFit(100);
  //toydata->Add(toydata1);

  //Generate and fit Ntoys samples of Poisson(nExpected) events, sig = nsig, bkg = nbkg
  //toydata->generateAndFit(Ntoys);

  //cout << toydata->fitParDataSet().get(10)->getRealValue("nsig_sb") << endl;
  //toydata->fitParDataSet().get(10)->Print("v");
  //toydata->fitParams(10)->Print("v");

  cout << "number of generated sig: " << nsig.getVal() << endl;
  cout << "number of fitted sig: " << nsig_sb.getVal() << " +/- " << nsig_sb.getError() << endl;
  //cout << "number of generated bkg: " << nbkg.getVal() << endl;

  //make plots
  //RooPlot* frame1 = toydata->plotParam(nsig, Bins(40)) ;
  //RooPlot* frame2 = toydata->plotError(nsig, Bins(40)) ;
  //RooPlot* frame3 = toydata->plotPull(nsig,Bins(40),FitGauss(kTRUE));

  
  //TH1F *hpull = new TH1F("hpull", "pull of nsig", 100, -10, 10);
  TH1F *hpull = new TH1F("hpull", "pull of nsig", 100, -0.5, 0.5);
  for (int i = 0; i < Ntoys; i++) {
    float pull_sig = 0.;
    RooRealVar* sig_result = (RooRealVar*) toydata->fitResult(i)->floatParsFinal().find(nsig_sb);
    //pull_sig = (sig_result->getValV() - nsig.getVal()) / sig_result->getError();
    pull_sig = (sig_result->getValV() - nsig.getVal()) /nsig.getVal();
    //RooRealVar* nkg_result = (RooRealVar*) toydata->fitResult(i)->floatParsFinal().find(nbkg);
    //double pull_sig = (sig_result->getValV() - nbkg.getVal()) / sig_result->getError();
    //cout << "fitted signal: " << sig_result->getValV() << "\t fitted bkg: " << nkg_result->getValV() << "\t pull: " << pull_sig << endl;
    hpull->Fill(pull_sig);
    //cout << "genpar datset: " <<  toydata->genData(i)->sumEntries() << endl;

  }

  
  gStyle->SetPalette(1) ;
  gStyle->SetOptStat(0) ;

  TCanvas* c1 = new TCanvas("nsig","nsig",600,500) ;
  c1->cd();
  //gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.4) ; frame1->Draw() ;

  TCanvas* c2 = new TCanvas("err","err",600,500) ;
  c2->cd();
  //gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.4) ; frame2->Draw() ;
  //gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0111);
  gStyle->SetStatY(0.95);
  gStyle->SetStatX(0.95);
  gStyle->SetStatW(0.18);
  gStyle->SetStatH(0.2);

  
  TCanvas* c3 = new TCanvas("pull","pull",600,500) ;
  c3->cd();
  //gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.4) ; frame3->Draw() ;
  hpull->GetYaxis()->SetTitle("Entries/0.5");
  hpull->GetXaxis()->SetTitle("pull of no. signal");
  hpull->Draw("e");
  TF1 *f1 = new TF1("gaus", "gaus", -3, 3);
  f1->SetLineColor(kBlue);
  f1->SetLineWidth(2);
  hpull->Fit(f1);
  

  TString dir = "plots/sibebandSel/SignalPull/";
  TString outname = "Signal_Pull_toyMC";

  if (barrel) outname += "_EB";
  else outname += "_EE";

  outname += Form("_PhoPt%dto%d", (int)minpt, (int)maxpt);

  outname += ".pdf";

  //c->SaveAs(dir + outname);


}
