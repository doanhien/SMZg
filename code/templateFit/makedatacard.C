#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCut.h"

#include <iostream>
#include <fstream>

void makedatacard(bool ele = true, bool barrel = true, TString cat = "pt", float minpt = 15, float maxpt = 20) {

  TString inf_da;

  //inf_da = "../ana_jet/fitting/minitrees/DoubleEG_DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_NJet_JetLooseID_Eta2p4.root";
  //inf_da = "../../ana_jet/fitting/minitrees/DoubleEG_DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_isPVGood.root";
  inf_da = "../../ana_jet/fitting/minitrees/DoubleEG_DoubleMu_Run2016Full_FebReminiAOD.root";

  cout << "input data file: " << inf_da.Data() << endl;

  TFile *fda = new TFile(inf_da, "read");
  //TFile *fsig_templ = new TFile("minitrees/Zg_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_allgenpho.root", "read");
  //TFile *fsig_templ = new TFile("../ana_jet/fitting/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr_NJet_JetLooseID_Eta2p4.root", "read");
  //TFile *fsig_templ = new TFile("../../ana_jet/fitting/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr_isPVGood_allSFs.root", "read");
  TFile *fsig_templ = new TFile("../../ana_jet/fitting/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr.root", "read");
  TFile *fbkg_templ = new TFile(inf_da, "read");

  TTree *tda = (TTree*)  fda->Get("outtree");
  TTree *tsig_templ = (TTree*) fsig_templ->Get("outtree");
  TTree *tbkg_templ = (TTree*) fbkg_templ->Get("outtree");

  TH1F *data_obs = new TH1F("data_obs", "hist from data", 10, -1, 1);
  TH1F *hsignal = new TH1F("hsignal", "signal template", 10, -1, 1);
  TH1F *background = new TH1F("background", "background template", 10, -1, 1);

  TH1F *signal_sigmaUp = new TH1F("signal_sigmaUp", "signal template minus 1sigma", 10, -1, 1);
  TH1F *signal_sigmaDown = new TH1F("signal_sigmaDown", "signal template plus 1sigma", 10, -1, 1);

  data_obs->Sumw2();
  hsignal->Sumw2();
  background->Sumw2();
  signal_sigmaUp->Sumw2();
  signal_sigmaDown->Sumw2();

  //TCut ptcut = Form("gamma_pt>%f && gamma_pt<%f", minpt, maxpt);
  TCut ptcut = Form("boss_pt>%f && boss_pt<%f", minpt, maxpt);
  TCut masscut = Form("boss_mass>%f && boss_mass<%f", minpt, maxpt);

  TCut cut = "z_charge==0 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4 && lept1_pt>20. && lept0_pt>25 && gamma_pt>=20";
  TCut cut_sb = "z_charge==0 && abs(lept0_eta)<2.4 && abs(lept1_eta)<2.4  && lept1_pt>20. && lept0_pt>25 && gamma_pt>=20";
  cut_sb += "(leptType==11 && trig_Ele23_Ele12==1) || (leptType==13 && trig_Mu17_Mu8==1 && pair_dPhi>70)";

  cut += "nselJet==0";
  cut_sb += "nselJet==0";

  if (ele) {
    cut += "leptType==11 && trig_Ele23_Ele12==1";
  }
  else {
    cut += "leptType==13 && trig_Mu17_Mu8==1 && pair_dPhi>70";
  }


  //side-band window from 24/12/2018
  if (barrel) {                                                                                                                                 
    cut += "isEB==1 && gamma_ChIso<2.";
    cut_sb += "isEB==1 && gamma_ChIso>7. && gamma_ChIso<13."; 
  }
  else {
    cut += "isEE==1 && gamma_ChIso<1.5";
    cut_sb += "isEE==1 && gamma_ChIso>6. && gamma_ChIso<14.";
  }

  if (cat == "pt") { 
    cut += ptcut;
    cut_sb += ptcut;
  }
  else {
    cut += masscut;
    cut_sb += masscut;
  }


  cout << "signal region: " << endl;
  cut.Print();
  cout << "side band region: " << endl;
  cut_sb.Print();

  tda->Draw("gamma_ssmva >> data_obs", cut, "goff");

  //cut *= "puweigj_65nb*genWeight";
  //if (ele) cut *= "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF";
  //else cut *= "puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  if (ele) cut *= "puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF";
  else cut *= "puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF";

  tsig_templ->Draw("gamma_ssmva >> hsignal", cut, "goff");

  tbkg_templ->Draw("gamma_ssmva >> background", (cut_sb), "goff");

  if (barrel) {
    tsig_templ->Draw("gamma_ssmva+0.06 >> signal_sigmaUp", cut, "goff");
    tsig_templ->Draw("gamma_ssmva-0.06 >> signal_sigmaDown", cut, "goff");
  
  }
  else {
    tsig_templ->Draw("gamma_ssmva+0.14 >> signal_sigmaUp", cut, "goff");
    tsig_templ->Draw("gamma_ssmva-0.14 >> signal_sigmaDown",cut, "goff");
  }

  TH1F *hgen_sig = new TH1F("hgen_sig", "hgen_sig", 10, -1, 1);
  TH1F *hgen_bkg = new TH1F("hgen_bkg", "hgen_bkg", 10, -1, 1);

  TH1F *hgen_sig_up = (TH1F*)signal_sigmaUp->Clone();
  TH1F *hgen_sig_down = (TH1F*)signal_sigmaDown->Clone();


  float bkg_bin9 = background->GetBinContent(9);
  float err_bkg_bin9 = background->GetBinError(9);
  //background->SetBinContent(9, bkg_bin9*2);
  //background->SetBinError(9, err_bkg_bin9*2);

  float sig_bin8 = hsignal->GetBinContent(7);
  float err_sig_bin8 = hsignal->GetBinError(7);

  //hsignal->SetBinContent(9, sig_bin8);
  //hsignal->SetBinError(9, err_sig_bin8);

  //hsignal->Scale(1./hsignal->Integral());
  //background->Scale(1./background->Integral());
  //signal_sigmaUp->Scale(1./signal_sigmaUp->Integral());
  //signal_sigmaDown->Scale(1./signal_sigmaDown->Integral());

  hntotweight = (TH1D*) fsig_templ->Get("hntotweight");
  float scale = 35900*47.34/hntotweight->Integral();

  hsignal->Scale(scale); //for bias
  cout << "entries from signal " << hsignal->Integral() << endl;

  ofstream fcard;
  TFile *fout;

  //TString outfname = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/Forbias/";
  TString outfname = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/exclusive/Forbias/";
  //TString outfname = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/signalTemplate_2sigma/";
  //TString outfname = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/Forbias/";
  TString fnamecard = outfname;

  if (cat == "mass") {
    outfname += "Mllg/datacard/data";
    fnamecard += "Mllg/datacard/datacard";
  }
  else {
    //outfname += "Ptg/datacard/data";
    //fnamecard += "Ptg/datacard/datacard";
    outfname += "Ptllg/datacard/data";
    fnamecard += "Ptllg/datacard/datacard";
  }
  

  TString sufix;

  if (cat == "pt" )
    //sufix= Form("_Pt%dto%d", (int) minpt, (int) maxpt);
    sufix= Form("_Ptllg%dto%d", (int) minpt, (int) maxpt);
  else
    sufix= Form("_Mllg%dto%d", (int) minpt, (int) maxpt);

  //if ( cat == "pt" ) {
    if (barrel) sufix += "_EB";
    else sufix += "_EE";
    //}
  if (ele) sufix += "_ele";
  else sufix += "_mu";

  outfname += sufix + ".root";
  fnamecard += sufix + ".txt";



  fout = new TFile(outfname, "recreate");
  fout->cd();

  data_obs->SetName("data_obs");
  data_obs->Write();
  hsignal->SetName("signal");
  hsignal->Write();
  signal_sigmaUp->Write();
  signal_sigmaDown->Write();
  background->Write();
  fout->Write();
  fout->Close();



  fcard.open(fnamecard);

  cout << "write data card" << endl;

  fcard << "#### pt:" << (int) minpt <<"-" << (int) maxpt << " GeV";
  if (barrel) fcard << " in EB,";
  else fcard << "EE,"; 
  if (ele) fcard << "ele channel" << endl;
  else fcard << "muon channel" << endl;
  //fcard  << "####" << endl;
  fcard << "------------------" << endl;
  fcard << "imax *" << endl;
  fcard << "jmax *" << endl;
  fcard << "kmax *" << endl;
  fcard << "------------------" << endl;
  //fcard << "shapes * * " << outfname.Data() << " $PROCESS" << endl;
  fcard << "shapes * * " << outfname.Data() << " $PROCESS" << " $PROCESS_$SYSTEMATIC" << endl;
  fcard << "------------------" << endl;
  fcard << "bin Zg" << endl;
  fcard << "observation " << data_obs->Integral() << endl;
  fcard << "----------------------------------------------" << endl;
  fcard << "bin \t \t Zg \t \t Zg" << endl;
  fcard << "process \t signal \t background" << endl;
  fcard << "process \t 0 \t \t 1" << endl;
  fcard << "rate \t \t " << hsignal->Integral() << "\t \t " << background->Integral() << endl;
  //fcard << "rate \t \t 1 \t \t 1" << endl;
  fcard << "----------------------------------------------" << endl;
  fcard << "sigma \t shape \t 1.0 \t-" << endl;
  //fcard << "SS rateParam Zg signal 1 [0," << data_obs->Integral() << "]" << endl;
  //fcard << "BB rateParam Zg background (" << data_obs->Integral() << "-@0) SS" << endl;
  //float bkg_exp = data_obs->Integral() - hsignal->Integral() + 1;
  fcard << "BB rateParam Zg background 1 [0," << data_obs->Integral() << "]" << endl;
  fcard.close();


  cout << "done !!!!!!!!!!!!!1" << endl;


}



