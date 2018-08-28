#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "untuplizer.h"
#include "tree.h"
#include "ElectronSelection.h"
#include "PhotonSelections.h"
//#include "cat.h"
#include "puweicalc.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>


void anaTnP_v2(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true)
{

  outputfile = "minitrees/";

  if (!mc)  outputfile = outputfile + "SingleEle";
  //else outputfile = outputfile + "data";

  if (inputfile.Contains("_Run2016B") ) outputfile += "_Run2016B";
  if (inputfile.Contains("_Run2016C") ) outputfile += "_Run2016C";
  if (inputfile.Contains("_Run2016D") ) outputfile += "_Run2016D";
  if (inputfile.Contains("_Run2016E") ) outputfile += "_Run2016E";
  if (inputfile.Contains("_Run2016F_ReMiniAOD1") ) outputfile += "_Run2016F_1";
  if (inputfile.Contains("_Run2016F_ReMiniAOD2") ) outputfile += "_Run2016F_2";
  if (inputfile.Contains("_Run2016G") ) outputfile += "_Run2016G";
  if (inputfile.Contains("_Run2016H_ReMiniAOD_03Feb_v2") ) outputfile += "_Run2016H_v2";
  if (inputfile.Contains("_Run2016H_ReMiniAOD_03Feb_v3") ) outputfile += "_Run2016H_v3";
  if (inputfile.Contains("_Zg_aMCatNLO/") ) outputfile += "_Zg_aMCatNLO";
  if (inputfile.Contains("_Zg_pt130/") ) outputfile += "_Zg_pt130";
  if (inputfile.Contains("DYJetsToLL_m50_MG") ) outputfile += "DYJets_MG";
  if (inputfile.Contains("DYJetsToLL_m50_aMCatNLO") ) outputfile += "DYJets_amcatnlo";


  outputfile += "_NoCalib_TnP_Zee_LooseSieie.root";

  TreeReader data(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *passingIdTree = new TTree("passingIdTree", "passing Id tree");

  inittree(passingIdTree);

  Long64_t nTotal(0) ;
  Long64_t nPassHLT(0) ;
  Long64_t nPassTag(0) ;
  Long64_t nFill(0) ;
  Long64_t ntot(0);

  PUWeightCalculator puCalcGJ;
  PUWeightCalculator puCalcGJ_69nb;
  PUWeightCalculator puCalcGJ_69p2nb;
  PUWeightCalculator puCalcGJ_65nb;
  PUWeightCalculator puCalcGJ_63nb;

  if ( data.HasMC() ) {
    //puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    puCalcGJ_69nb.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    puCalcGJ_69p2nb.Init("external/PU_histo_13TeV_GoldenJSON_69200nb.root");
    puCalcGJ_65nb.Init("external/PU_histo_13TeV_GoldenJSON_65550nb.root");
    puCalcGJ_63nb.Init("external/PU_histo_13TeV_GoldenJSON_63000nb_summer16.root");
  }

  Long64_t TotalEntries = data.GetEntriesFast();  
  for ( Long64_t ev = 0; ev < TotalEntries; ev++) {   
    //for ( Long64_t ev = 0; ev < 10000000; ev++) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);
    readggtree(data);

    nTotal++ ;

    if ( nEle < 2 ) continue ;

    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    isPVGood_ = isPVGood;
    run_ = run;
    event_ = event;
    lumi_ = lumis;

    //pfmet
    pfMet = pfMET;
    pfMetPhi = pfMETPhi;

    //PU weight for MC
    if (data.HasMC()) {
      float* puTrue = data.GetPtrFloat("puTrue");
      //puweigj = (float) puCalcGJ.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69nb = (float) puCalcGJ_69nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69p2nb = (float) puCalcGJ_69p2nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_65nb = (float) puCalcGJ_65nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_63nb = (float) puCalcGJ_63nb.GetWeight(run, puTrue[1]); // in-time PU

      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    //trigger
    //trig_Ele27_Eta2p1_WPTight = hlt>>1&1;
    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele23_Ele12 = hlt>>5&1;
    trig_DoubleEle33 = hlt>>10&1;

    //trig_Pho30 = hltPho>>1&1;
    trig_doublePho60 = hltPho >> 22&1;

    if (trig_doublePho60 == 1) nPassHLT++;

    vector<int> acc_tag_ele, acc_probe_ele;
    int nZ = 0;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 25.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      //if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data, i);
        if (!match) continue;
      }

      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1) acc_tag_ele.push_back(i);
      acc_probe_ele.push_back(i);
    }

    if (acc_tag_ele.size() < 1) continue;

    nTag = 0;

    for (vector<int>::iterator itag = acc_tag_ele.begin(); itag != acc_tag_ele.end(); itag++) {
      if (fabs(eleSCEta[*itag]) > 1.4442 && fabs(eleSCEta[*itag]) < 1.566) continue;
      if (fabs(eleSCEta[*itag]) > 2.1) continue;
      if ( elePt[*itag] < 30.) continue;

      //initial value;
      Tag_Pt = -99.;
      Tag_Eta = -99.;
      Tag_Phi = -99.;
      Tag_SCEta = -99.;
      TagIndex = -99;
      Tag_sieie = -99.;
      Tag_chIso = -99.;
      Tag_neuIso = -99.;
      Tag_phoIso = -99.;
      Tag_hoe = -99.;
      Tag_eop = -99.;
      Tag_ch = -99.;
      Tag_dEtaSeed = -99.;
      Tag_dEtaVtx = -99.;
      Tag_dPhiVtx = -99.;

      nProbe = 0;
      nPassId = 0;
      nFailId = 0;
      nPassIdIso = 0;
      nFailIdIso = 0;
      ProbeIndex = -99;
      Probe_Pt = -99.;
      Probe_Eta = -99.;
      Probe_Phi = -99.;
      Probe_SCEta = -99.;
      ProbeIndex = -99;
      Probe_sieie = -99.;
      Probe_chIso = -99.;
      Probe_neuIso = -99.;
      Probe_phoIso = -99.;
      Probe_hoe = -99.;
      Probe_eop = -99.;
      Probe_ch = -99.;
      Probe_dEtaSeed = -99.;
      Probe_dEtaVtx = -99.;
      Probe_dPhiVtx = -99.;
      Probe_miniIso = -99.;
      LooseID = -1;
      MediumID = -1;
      TightID = -1;
      LooseIDOnly = -1;
      Probe_PhoEleVeto = 0;
      Probe_PhoChIso = -99.;
      passPhoID_Zg = 0;

      Zm = -99.;
      Zeta = -99.;
      Zpt = -99.;

      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<float> ProbePt;
      vector<int> index;
      nProbe = 0;

      //for (vector<int>::iterator iele = acc_probe_ele.begin(); iele != acc_probe_ele.end(); iele++) {
      //if ( *itag == *iele) continue;
      //TLorentzVector vProbe ,vZ ;

      //vProbe.SetPtEtaPhiM(elePt[*iele], eleEta[*iele], elePhi[*iele],0.511*0.001);
      //vZ = vTag + vProbe;

      if (nPho < 1) continue;
      for (int ipho = 0; ipho < nPho; ipho++) {
	if (phoEt[ipho] < 15.) continue;
	if (fabs(phoSCEta[ipho]) > 2.5) continue;
	if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;

	//not match to tag ele
	if ( fabs(phoSCEta[ipho] - eleSCEta[*itag]) <= 0.05 || fabs(deltaPhi(phoSCPhi[ipho], eleSCPhi[*itag]) <= 0.05) ) continue;

	//match to gen ele
	if ( data.HasMC()) 
	  if (phoMatcher(data, ipho) != 1) continue;

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(phoEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	vZ = vTag + vProbe; 

	if (mc) {
	  if ( vZ.M() < 40 || vZ.M() > 140 ) continue ;  //for MC template
	} else {
	  if ( vZ.M() < 60 || vZ.M() > 120 ) continue ; //for data
	}
	  
	
	Zmass.push_back(vZ.M());
	ZEta.push_back(vZ.Eta());
	ZPt.push_back(vZ.Pt());
	//ProbePt.push_back(elePt[*iele]);
	//index.push_back(*iele);
	ProbePt.push_back(phoEt[ipho]);
	index.push_back(ipho);
      }
    
      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;

      //choose Z having closest to truth mass
      /*for (int i = 0; i < nZ; i++) {
	//cout << "Zmass: " << Zmass[i];
	float deltaM = fabs(91.19-Zmass[i]);
	if (minZ > deltaM) {
	  minZ = deltaM;
	  iProbe = index[i];
	  Zm = Zmass[i];
	  Zeta = ZEta[i];
	  Zpt = ZPt[i];
	}
	}*/

      //choose highest probe pt
      float maxpt = 0.;
      for (int i = 0; i < nZ; i++) {
	if (ProbePt[i] > maxpt) {
	  maxpt = ProbePt[i];
	  Zm = Zmass[i];
          Zeta = ZEta[i];
          Zpt = ZPt[i];
	  iProbe = index[i];
	}
      }

    
      if (iProbe < 0) continue;
      //if ( *itag == iProbe) continue;

      bool passId = false;
      bool passIdIso = false;

      nTag++;    
      Tag_Pt = elePt[*itag];
      Tag_Eta = eleEta[*itag];
      Tag_SCEta= eleSCEta[*itag];
      Tag_Phi = elePhi[*itag];
      Tag_sieie = eleSigmaIEtaIEta_Full5x5[*itag];
      Tag_chIso = elePFChIso[*itag];
      Tag_neuIso = elePFNeuIso[*itag];
      Tag_phoIso = elePFPhoIso[*itag];
      Tag_hoe = eleHoverE[*itag];
      Tag_eop = eleEoverP[*itag];
      Tag_dEtaSeed = eledEtaseedAtVtx[*itag];
      Tag_dEtaVtx = eledEtaAtVtx[*itag];
      Tag_dPhiVtx = eledPhiAtVtx[*itag];
      Tag_ch = eleCharge[*itag];
      Tag_IDMVA = eleIDMVA[*itag];
      ele1_MatchTrig27 = eleFiredSingleTrgs[*itag] >>12&1;
      mW_Tag_Met = sqrt(2*pfMet*elePt[*itag]*(1-cos(deltaPhi(pfMetPhi,elePhi[*itag]))));
      /*
      Probe_Pt = elePt[iProbe];
      Probe_Eta = eleEta[iProbe];
      Probe_SCEta = eleSCEta[iProbe];
      Probe_Phi = elePhi[iProbe];
      Probe_sieie = eleSigmaIEtaIEta_Full5x5[iProbe];
      Probe_chIso = elePFChIso[iProbe];
      Probe_neuIso = elePFNeuIso[iProbe];
      Probe_phoIso = elePFPhoIso[iProbe];
      Probe_hoe = eleHoverE[iProbe];
      Probe_eop = eleEoverP[iProbe];
      Probe_dEtaSeed = eledEtaseedAtVtx[iProbe];
      Probe_dEtaVtx = eledEtaAtVtx[iProbe];
      Probe_dPhiVtx = eledPhiAtVtx[iProbe];
      Probe_ch = eleCharge[iProbe];
      Probe_miniIso = elePFMiniIso[iProbe];

      LooseID = eleID[iProbe] >>1&1;
      MediumID = eleID[iProbe] >>2&1;
      TightID = eleID[iProbe] >>3&1;
      if ( pass_cutbased_80X(iProbe,0) ) LooseIDOnly = 1;
      else LooseIDOnly = 0;
      */

      Probe_Pt = phoEt[iProbe];
      Probe_Eta = phoEta[iProbe];
      Probe_SCEta = phoSCEta[iProbe];
      Probe_Phi = phoPhi[iProbe];
      Probe_PhoEleVeto = phoEleVeto[iProbe];
      //Probe_PhoChIso = phoPFChIso[iProbe];

      //rho-correction charged isolation 
      float chEA[] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
      int bin;
      if (phoSCEta[iProbe] < 1.0) bin = 0;
      else if (phoSCEta[iProbe] >= 1.0 && phoSCEta[iProbe] < 1.479) bin = 1;
      else if (phoSCEta[iProbe] >= 1.479 && phoSCEta[iProbe] < 2.0) bin = 2;
      else if (phoSCEta[iProbe] >= 2.0 && phoSCEta[iProbe] < 2.2) bin = 3;
      else if (phoSCEta[iProbe] >= 2.2 && phoSCEta[iProbe] < 2.3) bin = 4;
      else if (phoSCEta[iProbe] >= 2.3 && phoSCEta[iProbe] < 2.4) bin = 5;
      else bin = 6; // phoSCEta[iProbe] >= 2.4

      double corr = phoPFChIso[iProbe] - chEA[bin] * rho;
      double corrchIso = TMath::Max(corr, 0.);

      Probe_PhoChIso = corrchIso;

      if ( preselect_pho_Zg(iProbe) == 1 )
	if ( (fabs(phoSCEta[iProbe]) < 1.5 && phoPFChIso[iProbe]<2.) || (fabs(phoSCEta[iProbe]) > 1.5 && phoPFChIso[iProbe] <1.5) )
	  passPhoID_Zg = 1;
     

      passingIdTree->Fill();

    }

  }// end of loop
  cout << "total events = " << nTotal << endl ;
  cout << "pass HLT     = " << nPassHLT << endl ;
  cout << "pas tags     = " << nPassTag << endl ;
  cout << "filled       = " << nFill << endl ;
  cout << "endloop " << endl ;
  
  fo->Write();
  fo->Close();
  delete fo;

}
