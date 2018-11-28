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
//#include <TRatioPlot.h>

#include <iostream>
#include <fstream>

using namespace RooFit;
using namespace std;

void fitting (bool ele = true, bool barrel = true, TString cat = "pt", float minpt = 15., float maxpt = 20.) {

  TString inf_da;
  if (ele) inf_da = "../ana/minitrees/DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_LooseSieie.root";
  else inf_da = "../ana/minitrees/DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_LooseSieie.root";

  cout << "input data file: " << inf_da.Data() << endl;

  TFile *fda = new TFile(inf_da, "read");
  TFile *fsig_templ = new TFile("../ana/minitrees/Zg_aMCatNLO_Summer16_TMVA420_LooseSieie_checkgenmatch.root", "read");
  //TFile *fsig_templ = new TFile("../sideband_selection/output/bkgTemplate_summer16_420_LooseSieie_job_summer16_GJets_Pt_15To6000_TuneCUETP8M1_pythia8.root", "read");
  //TFile *fsig_templ = new TFile("../sideband_selection/output/bkgTemplate_summer16_420_LooseSieie_job_summer16_WGToLNuG_amcatnlo_ext3.root", "read");
  TFile *fbkg_templ = new TFile(inf_da, "read");


  TTree *tda = (TTree*)  fda->Get("outtree");
  TTree *tsig_templ = (TTree*) fsig_templ->Get("outtree");
  TTree *tbkg_templ = (TTree*) fbkg_templ->Get("outtree");

  /*
  TChain *tda_sb = new TChain("data");
  tda_sb->Add("../ana/minitrees/DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_LooseSieie.root/outtree");
  tda_sb->Add("../ana/minitrees/DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_LooseSieie.root/outtree");
  */

  TH1F *hda = new TH1F("hda", "hist from data", 20, -1, 1);
  TH1F *hsig_templ = new TH1F("hsig_templ", "signal template", 20, -1, 1);
  TH1F *hbkg_templ = new TH1F("hbkg_templ", "background template", 20, -1, 1);
  //TH1F *hbkg_templ = new TH1F("hbkg_templ", "background template", 18, -1, 0.8);


  TCut ptcut = Form("gamma_pt>%f && gamma_pt<%f", minpt, maxpt);
  TCut masscut = Form("boss_mass>%f && boss_mass<%f", minpt, maxpt);


  TCut cut = "z_charge==0";
  TCut cut_sb = "z_charge==0";

  if (ele) {
    cut += "leptType==11 && trig_Ele23_Ele12==1 && lept1_pt>20.";
    cut_sb += "leptType==11 && trig_Ele23_Ele12==1 && lept1_pt>20.";
  }
  else {
    cut += "leptType==13 && trig_Mu17_Mu8==1 && lept1_pt>20.";
    cut_sb += "leptType==13 && trig_Mu17_Mu8==1 && lept1_pt>20.";
  }

  if (barrel) {
    cut += "isEB==1 && gamma_ChIso<2.";
    cut_sb += "isEB==1 && gamma_ChIso>7. && gamma_ChIso<13.";
  }
  else {
    cut += "isEE==1 && gamma_ChIso<1.5";

    if ( minpt >= 25.) cut_sb += "isEE==1 && gamma_ChIso>5. && gamma_ChIso<11.";
    else cut_sb += "isEE==1 && gamma_ChIso>7. && gamma_ChIso<13.";
  }

  if ( cat == "pt") {cut += ptcut; cut_sb += ptcut;}
  else if ( cat == "mass") {cut += masscut; cut_sb += masscut;}


  TCut cut_tot_sb = "(leptType==11 && trig_Ele23_Ele12==1 && lept1_pt>20.) || (leptType==13 && trig_Mu17_Mu8==1 && lept1_pt>20.)";
  if (barrel) 
    cut_tot_sb += "isEB==1 && gamma_ChIso>7. && gamma_ChIso<13.";
  else cut_tot_sb += "isEE==1 && gamma_ChIso>7. && gamma_ChIso<13.";

  cout << "--- cut signal----:" ;
  cut.Print();
  cout << "--- cut background----:";
  cut_sb.Print();

  tda->Draw("gamma_ssmva >> hda", cut);

  TCut w = "puweigj_65nb*genWeight";

  TCut cut_sig = "";
  if (barrel) 
    cut_sig += "isEB==1 && gamma_ChIso<2.";
  else 
    cut_sig += "isEE==1 && gamma_ChIso<1.5";

  cut_sig += ptcut;  

  tsig_templ->Draw("gamma_ssmva >> hsig_templ", cut);
  //tsig_templ->Draw("gamma_ssmva >> hsig_templ", cut_sig);  //use Wg or gjet for template
  tbkg_templ->Draw("gamma_ssmva >> hbkg_templ", cut_sb);
  //tda_sb->Draw("gamma_ssmva >> hbkg_templ", cut_tot_sb);




  RooRealVar gamma_ssmva("gamma_ssmva", "Photon MVA", -1, 1);

  RooDataHist datah("datah", "# events in data", gamma_ssmva, Import(*hda));
  RooDataHist dahsig("dahsig", "# events in sigal", gamma_ssmva, Import(*hsig_templ));
  RooDataHist dahbkg("dahbkg", "# events in background", gamma_ssmva, Import(*hbkg_templ));


  //template pdf from histogram
  RooHistPdf sigpdf("sigpdf", "template pdf for signal", gamma_ssmva, dahsig);
  RooHistPdf bkgpdf("bkgpdf", "template pdf for background", gamma_ssmva, dahbkg);

  RooRealVar nsig("nsig", "#signal event", 400., 0., 200000.);
  RooRealVar nbkg("nbkg", "#background event", 200., 0., 200000.);

  //model to data
  RooAddPdf model("model", "sig+bkg", RooArgList(sigpdf, bkgpdf), RooArgList(nsig, nbkg));

  //extended likelihood fit
  model.fitTo(datah, Extended());

  RooPlot *frame1 = gamma_ssmva.frame(Title("fit data with signal and background template"));

  datah.plotOn(frame1);
  //model.fitTo(data);
  model.fitTo(datah, Extended());
  model.plotOn(frame1);

  //model.paramOn(frame1, Format("NEU", AutoPrecision(2)), Layout(0.4, 0.8, 0.9), LineWidth(0), LineColor(0));
  datah.plotOn(frame1, RooFit::Name("datah"), MarkerSize(1.0));
  model.plotOn(frame1, RooFit::Name("model"), LineColor(kBlue), LineWidth(2));
  //model.plotOn(frame1, RooFit::Name("sigpdf"), Components("sigpdf"),LineColor(kRed), LineWidth(3), LineStyle(kDashed));
  model.plotOn(frame1, RooFit::Name("sigpdf"), Components("sigpdf"),LineColor(kCyan-3), FillColor(kCyan-3), LineWidth(1), DrawOption("F"));
  model.plotOn(frame1, RooFit::Name("bkgpdf"), Components("bkgpdf"),LineColor(kRed-4),LineWidth(1),FillColor(kRed-4), FillStyle(3002),DrawOption("F"));
  datah.plotOn(frame1, RooFit::Name("datah"), MarkerSize(1.0));


  cout  << "----------- extracted signal--------------" << nsig.getVal() << endl;
  //frame1->getAttText()->SetTextFont(42);
  //frame1->getAttText()->SetTextSize(0.032);

  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  //pad1->SetLogy();
  pad1->SetFillColor(0);

  //frame1->SetLabelColor(0, "X");
  frame1->SetLabelSize(0, "X");
  frame1->SetTitleOffset(1.5, "Y");
  frame1->SetTitleOffset(1.3, "X");
  frame1->Draw("same");
  TLegend *lg = new TLegend(0.38, 0.4, 0.84, 0.64);
  lg->AddEntry(frame1->findObject("datah"),Form("data    %.1f events", hda->GetEntries()),"pe");
  lg->AddEntry(frame1->findObject("model"),Form("fitted   %.1f events", nsig.getVal() + nbkg.getVal()),"l");
  lg->AddEntry(frame1->findObject("sigpdf"),Form("signal  %.1f #pm %.1f events", nsig.getVal(), nsig.getError()),"f");
  lg->AddEntry(frame1->findObject("bkgpdf"),"Background","f");
  lg->Draw();
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.04);
  tx.SetTextFont(42);
  if (ele) tx.DrawLatex(0.44, 0.84, "Z#gamma#rightarrow ee#gamma");
  else tx.DrawLatex(0.4, 0.84, "Z#gamma #rightarrow #mu#mu#gamma");
  if (barrel) tx.DrawLatex(0.4, 0.78, "0 < |#eta^{#gamma}| < 1.5");
  else tx.DrawLatex(0.4,0.78, "1.5 < |#eta^{#gamma}| < 2.5");
  tx.DrawLatex(0.4,0.72, Form("%d < p_{T}^{#gamma} < %d", (int)minpt, (int)maxpt));

  RooDataHist *pdfDataHis = model.generateBinned(gamma_ssmva,0,true);
  TH1F *h_pdfDataHis = (TH1F*) pdfDataHis->createHistogram("h_pdfDataHis",gamma_ssmva, Binning(20));
  TH1F *h_fromRooData = (TH1F*)datah.createHistogram("h_fromRooData", gamma_ssmva, Binning(20));

  Double_t chi2=0.;
  //chi2=h_fromRooData->Chi2Test(h_pdfDataHis,"P CHI2");
  chi2=h_fromRooData->Chi2Test(h_pdfDataHis,"CHI2/NDF");
  //tx.DrawLatex(0.3, 0.8, Form("#chi^{2}/ndf = %.2f/20", frame1->chiSquare()));
  tx.DrawLatex(0.5, 0.35, Form("#chi^{2}/ndf = %.1f", frame1->chiSquare("model", "datah")));


  TH1F *hratio = (TH1F*) hda->Clone();
  //TH1F *hratio = (TH1F*) h_fromRooData->Clone();
  hratio->Divide(h_pdfDataHis);

  gStyle->SetOptStat(0);

  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  
  hratio->GetYaxis()->SetTitleOffset(1.6);
  hratio->SetMinimum(-6);
  hratio->SetTitleSize(0.14, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitle("data/ fit");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0., 2.);
  hratio->GetYaxis()->SetNdivisions(8);
  hratio->Draw();

  
  TLine *l = new TLine(-1,1,1,1);
  l->SetLineWidth(2);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");
  tx.SetTextFont(42);
  tx.SetTextSize(0.1);

  cout << "chi2 of fitting: " << frame1->chiSquare("model", "datah") << endl;
  cout << "chi2 of fitting: " << frame1->chiSquare() << endl;



  //file store signal yield

  ofstream fout;
  string filename;
  filename = "results/signal_yield";
  if (barrel) filename += "_EB";
  else filename += "_EE";
  if (ele) filename += "_ele.txt";
  else filename += "_mu.txt";

  //fout.open(filename, ios::app);

  //fout << (int) minpt << " < pt < " << (int) maxpt << "\t " << nsig.getVal() << " +/- " << nsig.getError() << endl;
  
  /*
  TCanvas *c2 = new TCanvas("c2", "c2", 650, 550);
  c2->cd();
  hbkg_templ->SetLineColor(2);
  hbkg_templ->Draw();
  hbkg_templ->Draw("same");
  cout << "bkg template: " << hbkg_templ->GetEntries() << endl;
  cout << "cut side band: " << endl;
  cut_sb.Print();
  */

  TString dir = "plots/fitting/private_code/";
  TString outname = "data_SigTempl_Zg_BkgTempl_dataSB";
  outname += Form("_PhoPt%dto%d", (int)minpt, (int)maxpt);
  if (barrel) outname += "_EB";
  else outname += "_EE";
  if (ele) outname += "_ele";
  else outname += "_mu";

  outname += ".pdf";

  //c1->SaveAs(dir + outname);

  //cout << "entries from data: " << hda->GetEntries() << endl;

}
