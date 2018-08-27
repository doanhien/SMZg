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

void xdataMCPlot(bool ele = true, TString sample = "data") {

  float lumi = 35900.; //whole 2016
  //float lumi = 5767; //2016B
  //float lumi = 2646; //2016C
  //float lumi = 4353; //2016D
  //float lumi = 3985; //2016E
  //float lumi = 3160; //2016F
  //float lumi = 7539; //2016G
  //float lumi = 8762; //2016H



  TString dir = "../ana/minitrees/";
  TString infname;
  if (sample.Contains("Zg") ) infname = "ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root";
  //if (sample.Contains("Zg") ) infname = "Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  //if (sample.Contains("Zg") ) infname = "Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root";
  //else if (sample.Contains("ZJet") ) infname = "ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_allfsr.root";
  else if (sample.Contains("ZJet") ) infname = "ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_genmatch_dr0p1.root";
  else if (sample.Contains("TT_") ) infname = "TT_Powheg_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("TTToSemilepton") ) infname = "TTToSemilepton_powheg_SortedPt.root";


  else if (sample.Contains("WWTo2L2Nu") ) infname = "WWTo2L2Nu_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("WWToLNuQQ") ) infname = "WWToLNuQQ_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("WZTo3LNu") ) infname = "WZTo3LNu_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("WZTo2L2Q") ) infname = "WZTo2L2Q_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("ZZTo2L2Nu") ) infname = "ZZTo2L2Nu_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("ZZTo2L2Q") ) infname = "ZZTo2L2Q_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  else if (sample.Contains("ZZTo4L") ) infname = "ZZTo4L_Summer16_TMVA420_UpTp6000_5VarCorr.root";

  else {
    if (ele) infname = "DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_nodoublecount.root";
    else infname = "DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_nodoublecount.root";
  }

  TFile *fin = new TFile(dir + infname, "READ");
  //TTree *tzg = (TTree*) fin->Get("tZ");
  TTree *tzg = (TTree*) fin->Get("outtree");

  cout << "processing sample: " << (dir+infname).Data() << endl;
  TString outfname;
  TString pref;
  if (ele) pref = "Ele_EB_";
  else pref = "Mu_EB_";
  TFile *fout = new TFile("histo/Kinematics/" + pref + infname , "recreate");

  //binning for lepton Pt
  float minPt, maxPt = 250.;
  float minPt_trail, maxPt_trail = 245;
  int binPt, binPt_trail;
  if (ele) {
    minPt = 25;
    minPt_trail = 20;
  }
  else {
    minPt = 20;
    maxPt = 245;
    minPt_trail = 20;
    maxPt_trail = 245;
  }

  binPt = (int) (maxPt-minPt)/3;
  binPt_trail = (int) (maxPt_trail-minPt_trail)/3;

  float minEta;
  if (ele) minEta = 2.5;
  else minEta = 2.4;

  const int nbinPt = 8;
  float phoPt[nbinPt+1] = {15, 20, 25, 30, 35, 45, 55, 65, 75};
  //define histogram
  TH1F *hlept1_Pt = new TH1F("hlept1_Pt", "; leading p_{T}^{e} (GeV/c); Entries/10 GeV", binPt, minPt, maxPt);
  TH1F *hlept2_Pt = new TH1F("hlept2_Pt", "; trailing p_{T}^{e} (GeV/c); Entries/4 GeV", binPt_trail, minPt_trail, maxPt_trail);
  TH1F *hlept1_Eta = new TH1F("hlept1_Eta", "; #eta^{e}; a.u", 20, -minEta, minEta);
  TH1F *hlept2_Eta = new TH1F("hlept2_Eta", "; #eta^{e}; a.u", 20, -minEta, minEta);
  TH1F *hZmass = new TH1F("hZmass", "; m_{e^{+}e^{-}} (GeV/c^{2}; Entries/2 GeV", 50, 50, 150);
  TH1F *hZmass_ss = new TH1F("hZmass_ss", "; m_{e^{+}e^{-}} (GeV/c^{2}; Entries/2 GeV", 50, 50, 150);
  //TH1F *hphoEt = new TH1F("hphoEt", "; p_{T}^{#gamma} (GeV/c); Entries/5 GeV", 45, 15, 105);
  TH1F *hphoEt = new TH1F("hphoEt", "; p_{T}^{#gamma} (GeV/c); Entries/5 GeV", nbinPt, phoPt);
  TH1F *hphoEta = new TH1F("hphoEta", "; #eta^{#gamma}; Entries/0.2", 30, -3, 3);
  TH1F *hmllg = new TH1F("hmllg", "; m_{X} (GeV/c^{2}); Events/2GeV", 65, 50, 180);
  TH1F *hdRPhoLep = new TH1F("hdRPhoLep", "; #deltaR(#gamma,l); Events", 50, 0, 5);
  TH1F *hnVtx = new TH1F("hnVtx", "", 50, 0, 50);
  TH1F *hssmva = new TH1F("hssmva", "", 20, -1, 1);

  //TCut cut = "z_mass>50 && z_charge==0 && ((fabs(gamma_eta)<1.5 && gamma_ChIso<2.) || (fabs(gamma_eta)>1.5 && gamma_ChIso<1.5))";
  TCut cut = "z_mass>50 && z_charge==0 && gamma_ChIso<2.";

  if (ele) cut += "leptType==11 && lept1_pt>20 && trig_Ele23_Ele12 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";
  else cut += "leptType==13 && trig_Mu17_Mu8 && lept1_pt>20 && lept0_pt>20 && lept1_pt>20";

  cut += "isEB";

  TCut cutss = "z_charge!=0 && fabs(gamma_eta)<2.5";
  if (ele) cutss += "leptType==11 && trig_Ele23_Ele12 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";
  else cutss += "leptType==13 && trig_Mu17_Mu8";
  TCut cutMllg = cut + "boss_mass<120";
  TCut cutpho = cut+ "gamma_pt < 45.";

  cut.Print();

  TCut weight;
  if (ele) weight = "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF";
  else weight = "puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF";
  //if (ele) weight = "puweigj_65nb*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF";
  //else weight = "puweigj_65nb*lept0_SelSF*lept1_SelSF*gamma_SF";

  cout << "weight for MC: " << endl;
  weight.Print();

  //read data to histogram
  if (sample == "data" ) {
    tzg->Draw("lept0_pt >> hlept1_Pt", cut);
    tzg->Draw("lept1_pt >> hlept2_Pt", cut);
    tzg->Draw("lept0_eta >> hlept1_Eta", cut);
    tzg->Draw("lept1_eta >> hlept2_Eta", cut);
    tzg->Draw("z_mass >> hZmass", cut);
    tzg->Draw("z_mass >> hZmass_ss", cutss);
    tzg->Draw("gamma_pt >> hphoEt", cutpho);
    tzg->Draw("gamma_eta >> hphoEta", cutpho);
    tzg->Draw("boss_mass >> hmllg", cutMllg);
    tzg->Draw("deltaR_PhoLept1 >> hdRPhoLep", cut);
    tzg->Draw("nVtx >> hnVtx", cut);
    tzg->Draw("gamma_ssmva >> hssmva", cut);

    //cout << "zmass before subtracting ss: " << hZmass->Integral() << endl;
    //hZmass->Add(hZmass_ss, -1);
    //cout << "zmass after subtracting ss: " << hZmass->Integral() << endl;

    cout << "selected Z in data: " << hZmass->Integral() << endl;
  }
  else {
    tzg->Draw("lept0_pt >> hlept1_Pt", cut*weight);
    tzg->Draw("lept1_pt >> hlept2_Pt", cut*weight);
    tzg->Draw("lept0_eta >> hlept1_Eta", cut*weight);
    tzg->Draw("lept1_eta >> hlept2_Eta", cut*weight);
    tzg->Draw("z_mass >> hZmass", cut*weight);
    tzg->Draw("z_mass >> hZmass_ss", cutss*weight);
    tzg->Draw("gamma_pt >> hphoEt", cutpho*weight);
    tzg->Draw("gamma_eta >> hphoEta", cutpho*weight);
    tzg->Draw("boss_mass >> hmllg", cutMllg*weight);
    tzg->Draw("deltaR_PhoLept1 >> hdRPhoLep", cut*weight);
    tzg->Draw("nVtx >> hnVtx", cut*weight);
    tzg->Draw("gamma_ssmva >> hssmva", cut*weight);

    //hZmass->Add(hZmass_ss, -1);
    //scale to luminosity
    TH1F *htotwei_Zg = (TH1F*) fin->Get("hntotweight");
    double totEvwei = htotwei_Zg->Integral();
    float xs;
    if (sample.Contains("Zg")) xs = 122.7;
    //if (sample.Contains("Zg")) xs = 47.34;
    else if (sample.Contains("ZJet")) xs = 5943.2;
    else if (sample.Contains("TT_")) xs = 87.31;
    else if (sample.Contains("TTToSemilepton")) xs = 320.1;
    else if (sample.Contains("WWTo2L2Nu")) xs = 12.178;
    else if (sample.Contains("WWToLNuQQ")) xs = 49.997;
    else if (sample.Contains("WZTo3LNu")) xs = 4.42965;
    else if (sample.Contains("WZTo2L2Q")) xs = 5.595;
    else if (sample.Contains("ZZTo2L2Nu")) xs = 0.564;
    else if (sample.Contains("ZZTo2L2Q")) xs = 3.22;
    else if (sample.Contains("ZZTo4L")) xs = 1.212;


    cout << "No of event weight: " << totEvwei << endl;
    float w_Zg = (lumi*xs/totEvwei);
    cout << "scale to lumi: " << w_Zg << endl;

    hlept1_Pt->Scale(w_Zg);
    hlept2_Pt->Scale(w_Zg);
    hlept1_Eta->Scale(w_Zg);
    hlept2_Eta->Scale(w_Zg);
    hZmass->Scale(w_Zg);
    hZmass_ss->Scale(w_Zg);
    hphoEt->Scale(w_Zg);
    hphoEta->Scale(w_Zg);
    hmllg->Scale(w_Zg);
    hdRPhoLep->Scale(w_Zg);
    hnVtx->Scale(w_Zg);
    hssmva->Scale(w_Zg);

    cout << "events from mc: " << hZmass->Integral() << endl;

  }

  //hZmass->Write();
  fout->Write();
  fout->Close();
}
