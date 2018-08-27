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

void xInputVarPlot(bool eb = true, TString sample = "data") {

  gStyle->SetOptStat(0);

  float lumi = 35900.;

  TString dir = "../ana/minitrees_Zee/";
  TString infname;
  //if ( sample.Contains("data") ) infname = "SingleEle_Run2016_TnP_Zee_chisocor_NewSSCor.root";
  if ( sample.Contains("data") ) infname = "SingleEle_Run2016_TnP_Zee_chisocor_NoGenMatch_NoSSCor.root";

  //else if (sample.Contains("DYJet") ) infname = "DYJets_amcatnlo_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr_BDTUpto6000.root";
  //else if (sample.Contains("TT") ) infname = "TTbar_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr_BDTUpto6000.root";
  //else if (sample.Contains("WJet") ) infname = "WJetsToLNu_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr_BDTUpto6000.root";

  //correction from Zee with PhoPresel
  else if (sample.Contains("DYJet") ) infname = "DYJets_amcatnlo_TnP_Zee_BDT_Upto6000_Corr_6Var_norho.root";
  else if (sample.Contains("TT") ) infname = "TTbar_TnP_Zee_BDT_Upto6000_Corr_6Var_norho.root";
  else if (sample.Contains("WJet") ) infname = "WJetsToLNu_TnP_Zee_BDT_Upto6000_Corr_6Var_norho.root";

  //else if (sample.Contains("DYJet") ) infname = "DYJets_amcatnlo_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root";
  //else if (sample.Contains("TT") ) infname = "TTbar_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root";
  //else if (sample.Contains("WJet") ) infname = "WJetsToLNu_TnP_Zee_chisocor_NewSSCor_rebin_PhoPresel_v2.root";

  TFile *fin = new TFile(dir + infname, "READ");
  TTree *tzg = (TTree*) fin->Get("passingIdTree");

  cout << "processing sample: " << (dir+infname).Data() << endl;
  TString outfname;
  TString pref;
  if (eb) pref = "EB_NoCorr_";
  else pref = "EE_NoCorr_";
  TFile *fout = new TFile("histo/inputVar_Trainning/" + pref + infname , "recreate");

  
  //----------photon Id ---------//
  TH1F *hphoHoE = new TH1F("hphoHoE", "; H/E; a.u", 50, 0, 0.1);
  TH1F *hphosieie = new TH1F("hphosieie", "; #sigma_{i#etai#eta}; a.u", 50, 0., 0.05);
  TH1F *hphochiso = new TH1F("hphochiso", "; charged iso^{#gamma}; a.u", 20, 0, 2);
  TH1F *hphowchiso = new TH1F("hphowchiso", "; worst charged iso^{#gamma}; a.u", 30, 0, 15);
  TH1F *hphossmva = new TH1F("hphossmva", "; SSMVA; a.u", 25, -1, 1);

  // --------input for shower shape---------//
  TH1F *hetawidth = new TH1F("hetawidth", "", 50, 0., 0.1);
  TH1F *hphiwidth = new TH1F("hphiwidth", "", 50, 0., 0.1);
  TH1F *hscRawE   = new TH1F("hscRawE", "", 60, 0, 600);
  TH1F *hs4 = new TH1F("hs4", "", 20, 0., 1.2);
  TH1F *hr9 = new TH1F("hr9", "", 20, 0., 1.2);
  TH1F *hsieip = new TH1F("hsieip", "", 50, -0.001, 0.001);
  TH1F *hphophi = new TH1F("hphophi", "", 30, -3.14, 3.14);
  TH1F *hphoeta = new TH1F("hphoeta", "", 25, -2.5, 2.5);
  TH1F *hrho = new TH1F("hrho", "", 50, 0, 50);


  TCut cut = "Zm>70 && Zm<110 && Probe_Pt>15 && passPhoID_Zg==1";
  //TCut cut = "Zm>70 && Zm<110 && Probe_Pt>15 && LooseID==1";
  if (eb) cut += "fabs(Probe_SCEta)<1.5";
  else cut += "fabs(Probe_SCEta)>1.5";


  cut.Print();

  TCut weight = "puweigj_65nb*genWeight*Tag_RecoSF*Tag_SelSF*Probe_phoSF";
  //if (ele) weight = "puweigj_65nb*genWeight*Tag_RecoSF*Tag_SelSF*Probe_phoSF";
  //else weight = "puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF";
  //weight = "puweigj*genWeight";

  cout << "weight for MC: " << endl;
  weight.Print();

  //read data to histogram
  if (sample.Contains("da")) {
    //tzg->Draw("Probe_HoverE >> hphoHoE", cut);
    tzg->Draw("Probe_sieie >> hphosieie", cut);
    tzg->Draw("Probe_sieip >> hsieip", cut);
    tzg->Draw("Probe_PhoChIso >> hphochiso", cut);
    tzg->Draw("Probe_ChWorstIso >> hphowchiso", cut);
    tzg->Draw("Probe_ssmva >> hphossmva", cut);
    tzg->Draw("Probe_SCEtaWidth >> hetawidth", cut);
    tzg->Draw("Probe_scphiwidth >> hphiwidth", cut);
    tzg->Draw("Probe_s4Full5x5 >> hs4", cut);
    tzg->Draw("Probe_R9 >> hr9", cut);
    tzg->Draw("Probe_SCRawE >> hscRawE", cut);
    tzg->Draw("rho >> hrho", cut);
    tzg->Draw("Probe_Phi >> hphophi", cut);
    tzg->Draw("Probe_SCEta >> hphoeta", cut);

  }
  else {
    //tzg->Draw("Probe_HoverE >> hphoHoE", cut*weight);
    tzg->Draw("Probe_sieie >> hphosieie", cut*weight);
    tzg->Draw("Probe_sieip >> hsieip", cut*weight);
    tzg->Draw("Probe_PhoChIso >> hphochiso", cut*weight);
    tzg->Draw("Probe_ChWorstIso >> hphowchiso", cut*weight);
    tzg->Draw("Probe_ssmva >> hphossmva", cut*weight);
    tzg->Draw("Probe_SCEtaWidth >> hetawidth", cut*weight);
    tzg->Draw("Probe_scphiwidth >> hphiwidth", cut*weight);
    tzg->Draw("Probe_s4Full5x5 >> hs4", cut*weight);
    tzg->Draw("Probe_R9 >> hr9", cut*weight);
    tzg->Draw("Probe_SCRawE >> hscRawE", cut*weight);
    tzg->Draw("rho >> hrho", cut*weight);
    tzg->Draw("Probe_Phi >> hphophi", cut*weight);
    tzg->Draw("Probe_SCEta >> hphoeta", cut*weight);
    
    //scale to luminosity
    /*
    TH1F *htotwei_Zg = (TH1F*) fin->Get("hntotweight");
    double totEvwei = htotwei_Zg->Integral();
    float xs;
    if (sample.Contains("Zg")) xs = 122.7;
    else if (sample.Contains("ZJet")) xs = 5943.2;
    else if (sample.Contains("TT")) xs = 87.31;
    else if (sample.Contains("WWTo2L2Nu")) xs = 12.178;
    else if (sample.Contains("WWToLNuQQ")) xs = 49.997;
    else if (sample.Contains("WZTo3LNu")) xs = 4.42965;
    else if (sample.Contains("WZTo2L2Q")) xs = 5.595;
    else if (sample.Contains("ZZTo2L2Nu")) xs = 0.564;
    else if (sample.Contains("ZZTo2L2Q")) xs = 3.22;
    else if (sample.Contains("ZZTo4L")) xs = 1.212;
    */

    //cout << "No of event weight: " << totEvwei << endl;
    float w_Zg = 1;
    if (sample.Contains("DYJet")) w_Zg = 2.60897;
    else if (sample.Contains("TT")) w_Zg = 0.0399477;
    else if (sample.Contains("WJet")) w_Zg = 35900*6.072e+04/160712265;

    cout << "scale to lumi: " << w_Zg << endl;

    //hphoHoE->Scale(w_Zg);
    hphosieie->Scale(w_Zg);
    hsieip->Scale(w_Zg);
    hphochiso->Scale(w_Zg);
    hphowchiso->Scale(w_Zg);
    hphossmva->Scale(w_Zg);
    hetawidth->Scale(w_Zg);
    hphiwidth->Scale(w_Zg);
    hs4->Scale(w_Zg);
    hr9->Scale(w_Zg);
    hrho->Scale(w_Zg);
    hphophi->Scale(w_Zg);
    hphoeta->Scale(w_Zg);
    hscRawE->Scale(w_Zg);

    cout << "events after scale: " << hphosieie->Integral() << endl;

  }

  fout->Write();
  fout->Close();

}
