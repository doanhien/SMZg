#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TF1.h>
#include <TMath.h>

#include "untuplizer.h"
#include "treeTnP.h"
#include "ElectronSelection.h"
//#include "MuonSelection.h"
#include "PhotonSelections.h"
#include "puweicalc.h"


void anaTnP(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true) {

  outputfile = "minitrees/";

  if (!mc) outputfile += "SingleEle";

  if (inputfile.Contains("_Run2016B") ) outputfile += "_Run2016B";
  if (inputfile.Contains("_Run2016C") ) outputfile += "_Run2016C";
  if (inputfile.Contains("_Run2016D") ) outputfile += "_Run2016D";
  if (inputfile.Contains("_Run2016E") ) outputfile += "_Run2016E";
  if (inputfile.Contains("_Run2016F_SepRereco1") ) outputfile += "_Run2016F_SepRereco1";
  if (inputfile.Contains("_Run2016F_SepRereco2") ) outputfile += "_Run2016F_SepRereco2";
  if (inputfile.Contains("_Run2016G") ) outputfile += "_Run2016G";
  if (inputfile.Contains("_Run2016H_PRv2") ) outputfile += "_Run2016H_PRv2";
  if (inputfile.Contains("_Run2016H") ) outputfile += "_Run2016H";

  if (inputfile.Contains("DYJet") ) outputfile += "DYJets";


  outputfile += "_TnP_Zee.root";

  TreeReader data1(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *outtree = new TTree("outtree", "output tree");
  inittree(outtree);
  TTree *outtreepho = new TTree("outtreepho", "output tree for pho");
  inittree(outtreepho);

  TH1D *hntotweight = new TH1D("hntotweight", "", 2, 1, 3);

  Long64_t ev1 = 0;
  Long64_t ev2 = -1;

  if (ev2 < 0) ev2 = data1.GetEntriesFast();
  if (ev2 > data1.GetEntriesFast()) ev2 = data1.GetEntriesFast();

  //PUWeightCalculator puCalcGJ;
  PUWeightCalculator puCalcGJ_69nb;
  PUWeightCalculator puCalcGJ_69p2nb;
  PUWeightCalculator puCalcGJ_65nb;
  PUWeightCalculator puCalcGJ_63nb;

  if ( data1.HasMC() ) {
    puCalcGJ_69nb.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    puCalcGJ_69p2nb.Init("external/PU_histo_13TeV_GoldenJSON_69200nb.root");
    puCalcGJ_65nb.Init("external/PU_histo_13TeV_GoldenJSON_65550nb.root");
    puCalcGJ_63nb.Init("external/PU_histo_13TeV_GoldenJSON_63000nb_summer16.root");
  }

  double ntot = 0;

  for (Long64_t ev = ev1; ev < ev2; ++ev) {
    //for (Long64_t ev = ev1; ev < 20000000; ++ev) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data1.GetEntriesFast());

    data1.GetEntry(ev);
    readggtree(data1);

    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    run_ = run;
    event_ = event;
    lumi_ = lumis;

    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele35_WPTight = hlt>>3&1;
    trig_Ele23_Ele12 = hlt>>5&1;
    trig_DoubleEle33 = hlt>>10&1;
    trig_doublePho70 = hltPho >> 22&1;

    //PU reweighting for MC
    if (data1.HasMC()) {
      float* puTrue = data1.GetPtrFloat("puTrue");
      puweigj_69nb = (float) puCalcGJ_69nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69p2nb = (float) puCalcGJ_69p2nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_65nb = (float) puCalcGJ_65nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_63nb = (float) puCalcGJ_63nb.GetWeight(run, puTrue[1]); // in-time PU
      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    if (trig_Ele27_WPTight != 1) continue;
    //if (trig_Ele35_WPTight != 1) continue;
    //if (nEle < 2) continue;

    vector<int> acc_tag_ele, acc_probe_ele;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 10.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data1, i);
        if (!match) continue;
      }

      //tight id for electron
      if ( pass_cutbased_80X(i,2) && (eleFiredSingleTrgs[i] >>12&1) ==1) acc_tag_ele.push_back(i); 
      //if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>39&1)==1) acc_tag_ele.push_back(i); //match to filter of HLT_Ele35
      if ( pass_cutbased_80X(i,0) ) acc_probe_ele.push_back(i);
      //if (eleIDMVA[i] > 0.2) acc_probe_ele.push_back(i);


    }

    if (acc_tag_ele.size() < 1) continue;
    if (acc_probe_ele.size() < 1) continue;

    int nZ = 0, nZ_matchPho = 0;
    ne = 0;
    npair = 0;
    npaireg = 0;

    for (unsigned int itag = 0 ; itag < acc_tag_ele.size(); itag++) {

      //initial values
      ProbeIndex = -99.;
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      tag_ele_MatchTrig27 = -99;

      probe_elePt = -99.;
      probe_eleEta = -99.;
      probe_eleSCEta = -99.;
      probe_elePhi = -99.;
      probe_ele_MatchTrig_Leg1_23 = -99;
      Zm = -99.;
      Zeta = -99.;
      //Zpt = -99.;

      if ( elePt[acc_tag_ele[itag]] < 30.) continue;
      nTag++;
      TagIndex = acc_tag_ele[itag];

      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[acc_tag_ele[itag]], eleEta[acc_tag_ele[itag]], elePhi[acc_tag_ele[itag]], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      nProbe = 0;

      for (unsigned int iele = 0; iele < acc_probe_ele.size(); iele++) {
	if ( acc_tag_ele[itag] >= acc_probe_ele[iele]) continue; //leading electron is tag
	if ( elePt[acc_probe_ele[iele]] < 10.) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(elePt[acc_probe_ele[iele]], eleEta[acc_probe_ele[iele]], elePhi[acc_probe_ele[iele]],0.511*0.001);
	vZ = vTag + vProbe;
	if ( vZ.M() < 60 || vZ.M() > 120 ) continue ;

	Zmass.push_back(vZ.M());
        ZEta.push_back(vZ.Eta());
        ZPt.push_back(vZ.Pt());
        index.push_back(acc_probe_ele[iele]);
      }

      //nZ++;
      //nProbe++;

      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
        float deltaM = fabs(91.19-Zmass[i]);
        if (minZ > deltaM) {
          minZ = deltaM;
          iProbe = index[i];
          Zm = Zmass[i];
          Zeta = ZEta[i];
          //Zpt = ZPt[i];
        }
      }

      if (iProbe < 0) continue;
      if ( acc_tag_ele[itag] == iProbe) continue;

      tag_elePt = elePt[acc_tag_ele[itag]];
      tag_eleEta = eleEta[acc_tag_ele[itag]];
      tag_eleSCEta= eleSCEta[acc_tag_ele[itag]];
      tag_elePhi = elePhi[acc_tag_ele[itag]];
      tag_ele_MatchTrig27 = eleFiredSingleTrgs[acc_tag_ele[itag]] >>12&1;
      //tag_ele_MatchTrig35 = eleFiredSingleTrgs[acc_tag_ele[itag]] >>39&1;
	
      ProbeIndex = iProbe;
      probe_elePt = elePt[iProbe];
      probe_eleEta = eleEta[iProbe];
      probe_eleSCEta = eleSCEta[iProbe];
      probe_elePhi = elePhi[iProbe];
      probe_ele_MatchTrig_Leg1_23 = eleFiredSingleTrgs[iProbe] >> 2&1;
      npair++;

      outtree->Fill();
    }


    //leading electron -trailing photon
    for (unsigned int itag = 0 ; itag < acc_tag_ele.size(); itag++) {

      //initial values
      ProbeIndex = -99.;
      tag_elePt = -99.;
      tag_eleEta = -99.;
      tag_eleSCEta = -99.;
      tag_elePhi = -99.;
      tag_ele_MatchTrig27 = -99;
      
      probe_phoPt = -99.;
      probe_phoEta = -99.;
      probe_phoSCEta = -99.;
      probe_phoPhi = -99.;
      probe_phoEleVeto = 0;
      probe_phoChIso = -99.;
      Zmass_ = -99.;
      Zeta_ = -99.;
      Zy_ = -99.;
   
      if ( elePt[acc_tag_ele[itag]] < 30.) continue;
      nTag++;
      TagIndex = acc_tag_ele[itag];
      
      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[acc_tag_ele[itag]], eleEta[acc_tag_ele[itag]], elePhi[acc_tag_ele[itag]], 0.511*0.001);
   
      vector<float> Zmass, ZEta, ZPt;
      vector<int> index;
      
      //for (int iele = 0; iele < nEle; iele++) {
      //if (acc_tag_ele[itag] >= iele) continue;
      //if (elePt[iele] < 15.) continue;
      //if (fabs(eleSCEta[iele]) > 2.5) continue;
      //if (fabs(eleSCEta[iele]) > 1.4442 && fabs(eleSCEta[iele]) < 1.566) continue;
	
	if (nPho < 1) continue;
	for (int ipho = 0; ipho < nPho; ipho++) {
	  if (phoCalibEt[ipho] < 15.) continue;
	  if (fabs(phoSCEta[ipho]) > 2.5) continue;
	  if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;
	  if ( preselect_pho_Zg(ipho) != 1) continue;
	  if ( fabs(phoSCEta[ipho]) < 1.5 && phoPFChIso[ipho] > 2.) continue;
	  if ( fabs(phoSCEta[ipho]) > 1.5 && phoPFChIso[ipho] > 1.5) continue;

	  //if ( phoIDMVA[ipho] < 0.2) continue;   
	  //if (phoEleVeto[ipho] !=1 ) continue;
	  //if ( fabs(phoSCEta[ipho] - eleSCEta[iele]) >= 0.05 || fabs(deltaPhi(phoSCPhi[ipho], eleSCPhi[iele]) >= 0.05) ) continue;
	  
	  TLorentzVector vProbe_matchPho ,vZ ;
	  vProbe_matchPho.SetPtEtaPhiM(phoCalibEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	  //vProbe_matchPho.SetPtEtaPhiM(phoEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	  vZ = vTag + vProbe_matchPho ;
	  if (vZ.M()< 60 || vZ.M() > 120 ) continue ;
	  
	  Zmass.push_back(vZ.M());
	  ZEta.push_back(vZ.Eta());
	  ZPt.push_back(vZ.Pt());
	  index.push_back(ipho);
	}
	//} //loop on ele
      
      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;
      for (int i = 0; i < nZ; i++) {
	float deltaM = fabs(91.19-Zmass[i]);
	if (minZ > deltaM) {
	  minZ = deltaM;
	  iProbe = index[i];
	  Zmass_ = Zmass[i];
	  Zeta_ = ZEta[i];
	  //Zy_matchPho = ZRapidity[i];
	}
      }
      
      if (iProbe < 0) continue;
      //if ( acc_tag_ele[itag] == iProbe) continue;
      //not select photon matched to leading electron
      if ( fabs(phoSCEta[iProbe] - eleSCEta[acc_tag_ele[itag]]) <= 0.05 || fabs(deltaPhi(phoSCPhi[iProbe], eleSCPhi[acc_tag_ele[itag]]) <= 0.05) ) continue;

      tag_elePt = elePt[acc_tag_ele[itag]];
      tag_eleEta = eleEta[acc_tag_ele[itag]];
      tag_eleSCEta= eleSCEta[acc_tag_ele[itag]];
      tag_elePhi = elePhi[acc_tag_ele[itag]];
      tag_ele_MatchTrig27 = eleFiredSingleTrgs[acc_tag_ele[itag]] >>12&1;
      
      probe_phoPt = phoEt[iProbe];
      probe_phoEta = phoEta[iProbe];
      probe_phoSCEta = phoSCEta[iProbe];
      probe_phoPhi = phoPhi[iProbe];
      probe_phoEleVeto = phoEleVeto[iProbe];
      probe_phoChIso = phoPFChIso[iProbe];

      //Zm_matchPho = vZ.M();
      //Zeta_matchPho = vZ.Eta();
      //Zy_matchPho = vZ.Rapidity();
      npaireg++;
      
      outtreepho->Fill();
    }
  }

  fo->Write();
  fo->Close();
  
}

