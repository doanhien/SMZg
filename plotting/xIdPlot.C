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

void xIdPlot(bool ele = true, TString sample = "data") {

  gStyle->SetOptStat(0);

  float lumi = 35900.;

  TString dir = "../ana/minitrees/";
  TString infname;
  if (sample.Contains("Zg") ) infname = "Zg_aMCatNLO_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("ZJet") ) infname = "ZJets_aMCatNLO_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("TT") ) infname = "TT_Powheg_Summer16_TMVA420_NoPhoSel.root";

  else if (sample.Contains("WWTo2L2Nu") ) infname = "WWTo2L2Nu_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("WWToLNuQQ") ) infname = "WWToLNuQQ_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("WZTo3LNu") ) infname = "WZTo3LNu_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("WZTo2L2Q") ) infname = "WZTo2L2Q_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("ZZTo2L2Nu") ) infname = "ZZTo2L2Nu_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("ZZTo2L2Q") ) infname = "ZZTo2L2Q_Summer16_TMVA420_NoPhoSel.root";
  else if (sample.Contains("ZZTo4L") ) infname = "ZZTo4L_Summer16_TMVA420_NoPhoSel.root";

  else {
    if (ele) infname = "DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_NoPhoSel.root";
    else infname = "DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_NoPhoSel.root";
  }

  TFile *fin = new TFile(dir + infname, "READ");
  TTree *tzg = (TTree*) fin->Get("outtree");

  cout << "processing sample: " << (dir+infname).Data() << endl;
  TString outfname;
  TString pref;
  if (ele) pref = "Ele_EE_";
  else pref = "Mu_EE_";
  TFile *fout = new TFile("histo/Id/" + pref + infname , "recreate");

  //data
  //---lepton Id------//
  TH1F *hlepd0_lead = new TH1F("hlepd0_lead", "; D0; a.u", 100, -0.1, 0.1);
  TH1F *hlepdz_lead = new TH1F("hlepdz_lead", "; DZ; a.u", 100, -0.4, 0.4);
  TH1F *hSIP_lead = new TH1F("hSIP_lead", "; SIP; a.u", 20, 0, 5);
  TH1F *hchiso_lead = new TH1F("hchiso_lead", "; charged iso; a.u", 100, 0, 10);
  TH1F *hphoiso_lead = new TH1F("hphoiso_lead", "; photon iso; a.u", 100, 0, 10);
  TH1F *hneuiso_lead = new TH1F("hneuiso_lead", "; neutral iso; a.u", 100, 0, 10);
  TH1F *hmva_lead = new TH1F("hmva_lead", "; MVA; a.u", 100, -1., 1.);
  TH1F *hsigEOverE_lead = new TH1F("hsigEOverE_lead", "; TrkPtError/TrkPt; a.u", 30, 0, 0.3);
  TH1F *hmuSta = new TH1F("hmuSta", "; muon stations; a.u", 8, 0, 8);
  TH1F *hmuPixhit = new TH1F("hmuPixhit", "; muon stations; a.u", 15, 0, 15);
  TH1F *hmutrkLayer = new TH1F("hmutrkLayer", "; muon stations; a.u", 20, 0, 20);

  TH1F *hlepd0_trail = new TH1F("hlepd0_trail", "; D0; a.u", 100, -0.1, 0.1);
  TH1F *hlepdz_trail = new TH1F("hlepdz_trail", "; DZ; a.u", 100, -0.4, 0.4);
  TH1F *hSIP_trail = new TH1F("hSIP_trail", "; SIP; a.u", 20, 0, 5);
  TH1F *hchiso_trail = new TH1F("hchiso_trail", "; charged iso; a.u", 100, 0, 10);
  TH1F *hphoiso_trail = new TH1F("hphoiso_trail", "; photon iso; a.u", 100, 0, 10);
  TH1F *hneuiso_trail = new TH1F("hneuiso_trail", "; neutral iso; a.u", 100, 0, 10);
  TH1F *hmva_trail = new TH1F("hmva_trail", "; MVA; a.u", 100, -1., 1.);
  TH1F *hsigEOverE_trail = new TH1F("hsigEOverE_trail", "; TrkPtError/TrkPt; a.u", 30, 0, 0.3);
  
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
  TH1F *hrho = new TH1F("hrho", "", 50, 0, 50);



  //TCut cut = "z_mass>50 && z_charge==0 && gamma_pt>15. && fabs(gamma_eta)<2.5 && ((fabs(gamma_eta)<1.5 && gamma_ChIso<2.) || (fabs(gamma_eta)>1.5 && gamma_ChIso<1.5))";
  TCut cut = "z_mass>50 && z_charge==0 && gamma_pt>15. && fabs(gamma_sceta)>1.5 && gamma_ChIso<15.0";
  if (ele) cut += "leptType==11 && lept0_pt>25. && lept1_pt>20. && trig_Ele23_Ele12";
  else cut += "leptType==13 && lept0_pt>20. && lept1_pt>20. && trig_Mu17_Mu8";

  //if (barrel) cut += "fabs(gamma_eta)<2.5";
  //else cut += "fabs(gamma_eta)>1.5";

  cut.Print();

  TCut weight;
  if (ele) weight = "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF";
  else weight = "puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF";
  //else  weight = "puweigj*genWeight";

  cout << "weight for MC: " << endl;
  weight.Print();

  //read data to histogram
  if (sample.Contains("da")) {
    tzg->Draw("lept0_d0 >> hlepd0_lead", cut);
    tzg->Draw("lept0_dz >> hlepdz_lead", cut);
    tzg->Draw("lept0_SIP >> hSIP_lead", cut);
    tzg->Draw("lept0_chiso >> hchiso_lead", cut);
    tzg->Draw("lept0_phoiso >> hphoiso_lead", cut);
    tzg->Draw("lept0_neuiso >> hneuiso_lead", cut);
    tzg->Draw("lept0_mva >> hmva_lead", cut);
    tzg->Draw("lept0_sigEOverE >> hsigEOverE_lead", cut);
    tzg->Draw("lept0_nMatchedStations >> hmuSta", cut);
    tzg->Draw("lept0_nValidPixelHits >> hmuPixhit", cut);
    tzg->Draw("lept0_nTrackerLayers >> hmutrkLayer", cut);
    
    tzg->Draw("lept1_d0 >> hlepd0_trail", cut);
    tzg->Draw("lept1_dz >> hlepdz_trail", cut);
    tzg->Draw("lept1_SIP >> hSIP_trail", cut);
    tzg->Draw("lept1_chiso >> hchiso_trail", cut);
    tzg->Draw("lept1_phoiso >> hphoiso_trail", cut);
    tzg->Draw("lept1_neuiso >> hneuiso_trail", cut);
    tzg->Draw("lept1_mva >> hmva_trail", cut);
    tzg->Draw("lept1_sigEOverE >> hsigEOverE_trail", cut);
    
    tzg->Draw("gamma_HoverE >> hphoHoE", cut);
    tzg->Draw("gamma_sigmaIetaIeta >> hphosieie", cut);
    tzg->Draw("gamma_ChIso >> hphochiso", cut);
    tzg->Draw("gamma_ChWorstIso >> hphowchiso", cut);
    tzg->Draw("gamma_ssmva >> hphossmva", cut);
    tzg->Draw("gamma_SCEtaWidth >> hetawidth", cut);
    tzg->Draw("gamma_scphiwidth >> hphiwidth", cut);
    tzg->Draw("gamma_s4Full5x5 >> hs4", cut);
    tzg->Draw("gamma_R9 >> hr9", cut);
    tzg->Draw("gamma_sigmaIetaIphi >> hsieip", cut);
    tzg->Draw("gamma_SCRawE >> hscRawE", cut);
    tzg->Draw("rho >> hrho", cut);
    tzg->Draw("gamma_phi >> hphophi", cut);

  }
  else {
    tzg->Draw("lept0_d0 >> hlepd0_lead", cut*weight);
    tzg->Draw("lept0_dz >> hlepdz_lead", cut*weight);
    tzg->Draw("lept0_SIP >> hSIP_lead", cut*weight);
    tzg->Draw("lept0_chiso >> hchiso_lead", cut*weight);
    tzg->Draw("lept0_phoiso >> hphoiso_lead", cut*weight);
    tzg->Draw("lept0_neuiso >> hneuiso_lead", cut*weight);
    tzg->Draw("lept0_mva >> hmva_lead", cut*weight);
    tzg->Draw("lept0_sigEOverE >> hsigEOverE_lead", cut*weight);
    tzg->Draw("lept0_nMatchedStations >> hmuSta", cut*weight);
    tzg->Draw("lept0_nValidPixelHits >> hmuPixhit", cut*weight);
    tzg->Draw("lept0_nTrackerLayers >> hmutrkLayer", cut*weight);
    
    tzg->Draw("lept1_d0 >> hlepd0_trail", cut*weight);
    tzg->Draw("lept1_dz >> hlepdz_trail", cut*weight);
    tzg->Draw("lept1_SIP >> hSIP_trail", cut*weight);
    tzg->Draw("lept1_chiso >> hchiso_trail", cut*weight);
    tzg->Draw("lept1_phoiso >> hphoiso_trail", cut*weight);
    tzg->Draw("lept1_neuiso >> hneuiso_trail", cut*weight);
    tzg->Draw("lept1_mva >> hmva_trail", cut*weight);
    tzg->Draw("lept1_sigEOverE >> hsigEOverE_trail", cut*weight);
    
    tzg->Draw("gamma_HoverE >> hphoHoE", cut*weight);
    tzg->Draw("gamma_sigmaIetaIeta_rw >> hphosieie", cut*weight);
    tzg->Draw("gamma_ChIso >> hphochiso", cut*weight);
    tzg->Draw("gamma_ChWorstIso >> hphowchiso", cut*weight);
    tzg->Draw("gamma_ssmva >> hphossmva", cut*weight);
    tzg->Draw("gamma_SCEtaWidth_rw >> hetawidth", cut*weight);
    tzg->Draw("gamma_scphiwidth >> hphiwidth", cut*weight);
    tzg->Draw("gamma_s4Full5x5_rw >> hs4", cut*weight);
    tzg->Draw("gamma_R9_rw >> hr9", cut*weight);
    tzg->Draw("gamma_sigmaIetaIphi >> hsieip", cut*weight);
    tzg->Draw("gamma_SCRawE >> hscRawE", cut*weight);
    tzg->Draw("rho >> hrho", cut*weight);
    tzg->Draw("gamma_phi >> hphophi", cut*weight);

    //scale to luminosity
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


    cout << "No of event weight: " << totEvwei << endl;
    float w_Zg = (lumi*xs/totEvwei);
    cout << "scale to lumi: " << w_Zg << endl;

    hlepd0_lead->Scale(w_Zg);
    hlepdz_lead->Scale(w_Zg);
    hSIP_lead->Scale(w_Zg);
    hchiso_lead->Scale(w_Zg);
    hphoiso_lead->Scale(w_Zg);
    hneuiso_lead->Scale(w_Zg);
    hmva_lead->Scale(w_Zg);
    hsigEOverE_lead->Scale(w_Zg);
    hmuSta->Scale(w_Zg);
    hmuPixhit->Scale(w_Zg);
    hmutrkLayer->Scale(w_Zg);
    
    hlepd0_trail->Scale(w_Zg);
    hlepdz_trail->Scale(w_Zg);
    hSIP_trail->Scale(w_Zg);
    hchiso_trail->Scale(w_Zg);
    hphoiso_trail->Scale(w_Zg);
    hneuiso_trail->Scale(w_Zg);
    hmva_trail->Scale(w_Zg);
    hsigEOverE_trail->Scale(w_Zg);

    hphoHoE->Scale(w_Zg);
    hphosieie->Scale(w_Zg);
    hphochiso->Scale(w_Zg);
    hphowchiso->Scale(w_Zg);
    hphossmva->Scale(w_Zg);
    hetawidth->Scale(w_Zg);
    hphiwidth->Scale(w_Zg);
    hs4->Scale(w_Zg);
    hr9->Scale(w_Zg);
    hrho->Scale(w_Zg);
    hsieip->Scale(w_Zg);
    hphophi->Scale(w_Zg);
    hscRawE->Scale(w_Zg);

    cout << "events after scale: " << hphoHoE->Integral() << endl;

  }

  fout->Write();
  fout->Close();

}
