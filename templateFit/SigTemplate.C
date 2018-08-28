#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "../ana/puweicalc.h"
#include "../ana/untuplizer.h"
#include "PhotonSelection.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"


using namespace std;

void SigTemplate(const char** inpaths, int npaths, TString outpath = "minitree.root", Float_t xs=1, Float_t lumi=1, Int_t aMCatNLO=0) {

  /// inpaths  = array of paths to files with ggNtuplizer's TTrees;
  /// npaths   = size of inpaths;
  /// outpath  = path to output root file;
  /// xs       = cross section
  /// lumi     = luminosity
  /// aMCatNLO = 0 : not aMC@NLO, 1 : aMC@NLO (with negative weight) 

  // open tree(s) with events to be processed
  TreeReader data(inpaths, npaths);

  Float_t puwei_  = 1.; 
  Float_t mcwei_  = 1.;
  Float_t genwei_ = 1.;
  Float_t gamma_pt_, gamma_ssmva_, gamma_ChIso_;
  Float_t gamma_eta_, gamma_sceta_, gamma_phi_;
  Float_t gamma_PhoIso_, gamma_sigIetaIeta_;
  Float_t gamma_ChWorstIso_, gamma_HoverE_;
  Int_t isEB_, isEE_;
  Int_t isEle_, isMu_;
  Int_t leptType_;
  Float_t z_mass_;
  Int_t   mcMomPID_; 
  Float_t isogen_;
  Float_t boss_mass_;
  Float_t boss_pt_;

  Int_t nseljet;
  Float_t jetPt_[10];
  Float_t jetEta_[10];
  Float_t jetPhi_[10];
  
  // prepare output tree
  TFile* fo = TFile::Open(outpath.Data(), "RECREATE");
  if (!fo || fo->IsZombie())  FATAL("TFile::Open() failed");
  //fo->mkdir("ggNtuplizer")->cd();
  fo->cd();
  TTree *outtree = new TTree("outtree", "mini tree for building signal template");
  outtree->Branch("puwei",         &puwei_,       "puwei/F");
  outtree->Branch("mcwei",         &mcwei_,       "mcwei/F");
  outtree->Branch("genwei",        &genwei_,      "genwei/F");
  outtree->Branch("gamma_pt",      &gamma_pt_,    "gamma_pt/F"); 
  outtree->Branch("gamma_eta",     &gamma_eta_,   "gamma_eta/F"); 
  outtree->Branch("gamma_sceta",   &gamma_sceta_,   "gamma_sceta/F"); 
  outtree->Branch("gamma_phi",     &gamma_phi_,   "gamma_phi/F"); 
  outtree->Branch("gamma_ssmva",   &gamma_ssmva_, "gamma_ssmva/F");  
  outtree->Branch("gamma_ChIso",   &gamma_ChIso_, "gamma_ChIso/F");
  outtree->Branch("gamma_ChWorstIso", &gamma_ChWorstIso_, "gamma_ChWorstIso/F");
  outtree->Branch("gamma_PhoIso",  &gamma_PhoIso_, "gamma_PhoIso/F");
  outtree->Branch("gamma_sigIetaIeta", &gamma_sigIetaIeta_, "gamma_sigIetaIeta/F");
  outtree->Branch("gamma_HoverE",  &gamma_HoverE_,  "gamma_HoverE/F");

  outtree->Branch("isEB",        &isEB_,        "isEB/I"); 
  outtree->Branch("isEE",        &isEE_,        "isEE/I"); 
  outtree->Branch("isEle",       &isEle_,       "isEle/I"); 
  outtree->Branch("isMu",        &isMu_,        "isMu/I"); 
  outtree->Branch("leptType",    &leptType_,    "leptType/I");
  outtree->Branch("z_mass",      &z_mass_,      "z_mass/F");
  outtree->Branch("mcMomPID",    &mcMomPID_,    "mcMomPID/I");
  outtree->Branch("isogen",      &isogen_,      "isogen/F"); 
  outtree->Branch("boss_mass",   &boss_mass_,   "boss_mass/F");
  outtree->Branch("boss_pt",     &boss_pt_,     "boss_pt/F");

  outtree->Branch("nseljet",      &nseljet);
  outtree->Branch("jetPt",        &jetPt_,     "jetPt[nseljet]/F");
  outtree->Branch("jetEta",       &jetEta_,   "jetEta[nseljet]/F");
  outtree->Branch("jetPhi",       &jetPhi_,   "jetPhi[nseljet]/F");


  // pileup reweighting for MC
  PUWeightCalculator puCalc;
  if (data.HasMC()) puCalc.Init("../ana/external/PU_histo_13TeV_GoldenJSON_69200nb.root");

  // compute mc weights
  if (aMCatNLO == 1) {
    Float_t totalEvents = 0;
    for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
      
      data.GetEntry(ev);
      if (data.HasMC()) {
	float genWeight = data.GetFloat("genWeight");
	if (genWeight > 0) totalEvents++;
	else totalEvents--;
      }
    }
    mcwei_ = (totalEvents != 0) ? xs * lumi / totalEvents : 1;    
    cout<<"total events : "<<totalEvents<<endl;
  } else {
    cout<<xs<<" "<<lumi<<" "<<data.GetEntriesFast()<<endl;
    mcwei_ = xs * lumi / data.GetEntriesFast();
  }

  Int_t counts = 0;
  Int_t countNoPho = 0;  

  //for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
  for (Long64_t ev = 0; ev < 40000000; ev++) {
    if ((ev % 50000 ) == 0) cout << "processing event " << ev+1 << " in " << data.GetEntriesFast() << endl;
    data.GetEntry(ev);
    
    Int_t run_  = data.GetInt("run");
    Float_t rho = data.GetFloat("rho");
    Int_t lumis_ = data.GetInt("lumis");
    Long64_t event_ = data.GetLong64("event");
    //nVtx_  = data.GetInt("nVtx");
    
    //cout << "processing event: " << event_ << " of lumis: " << lumis_ << endl;
    // PU reweighting
    if (data.HasMC()) {
      float* puTrue = data.GetPtrFloat("puTrue");
      puwei_ = (float) puCalc.GetWeight(run_, puTrue[1]); // in-time PU
      
      float generatorWeight = data.GetFloat("genWeight");
      if (aMCatNLO == 1) genwei_ = (generatorWeight > 0) ? 1. : -1.;
      else               genwei_ = generatorWeight;
    }
    
    Int_t  nMC      = 0;
    Int_t* mcPID    = NULL;
    Int_t* mcMomPID = NULL;
    Float_t* mcPt     = NULL;
    Float_t* mcEta    = NULL;
    Float_t* mcPhi    = NULL;
    UShort_t* mcStatusFlag = NULL;
    Float_t* mcCalIsoDR03 = NULL;
    if (data.HasMC()) {
      nMC          = data.GetInt("nMC");
      mcPID        = data.GetPtrInt("mcPID");
      mcMomPID     = data.GetPtrInt("mcMomPID");
      mcPt         = data.GetPtrFloat("mcPt");
      mcEta        = data.GetPtrFloat("mcEta");
      mcPhi        = data.GetPtrFloat("mcPhi");
      mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
      mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");
    }

    TLorentzVector mcpho, pho;
    
    Int_t  nPho       = data.GetInt("nPho");
    Float_t* phoEt      = data.GetPtrFloat("phoCalibEt");
    Float_t* phoSCEta   = data.GetPtrFloat("phoSCEta");
    Float_t* phoEta     = data.GetPtrFloat("phoEta");
    Float_t* phoPhi     = data.GetPtrFloat("phoPhi");
    Float_t* phoPFChIso = data.GetPtrFloat("phoPFChIso");   
    Float_t* phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
    Float_t* phoPFChWorstIso = data.GetPtrFloat("phoPFChWorstIso");
    Float_t* phoHoverE = data.GetPtrFloat("phoHoverE");
    Float_t* phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");

    Int_t nEle        = data.GetInt("nEle");
    Float_t* elePt    = data.GetPtrFloat("elePt");
    Float_t* eleEta   = data.GetPtrFloat("eleEta");
    Float_t* eleSCEta = data.GetPtrFloat("eleSCEta");
    Float_t* elePhi   = data.GetPtrFloat("elePhi");
    Short_t* eleID    = data.GetPtrShort("eleIDbit");

    Int_t nMu         = data.GetInt("nMu");
    Float_t* muPt     = data.GetPtrFloat("muPt");
    Float_t* muEta    = data.GetPtrFloat("muEta");
    Float_t* muPhi    = data.GetPtrFloat("muPhi");
    //Short_t* muIDbit = data.GetPtrShort("muIDbit");
    Float_t* muPFChIso = data.GetPtrFloat("muPFChIso");
    Float_t* muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
    Float_t* muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
    Float_t* muPFPUIso = data.GetPtrFloat("muPFPUIso");

    ULong64_t HLTJet = data.GetLong64("HLTJet");
    Int_t nJet = data.GetInt("nJet");
    Float_t* jetPt = data.GetPtrFloat("jetPt");
    Float_t* jetEta = data.GetPtrFloat("jetEta");
    Float_t* jetPhi = data.GetPtrFloat("jetPhi");
    Float_t* jetEn = data.GetPtrFloat("jetEn");
    Float_t* jetRawPt = data.GetPtrFloat("jetRawPt");
    Float_t* jetRawEn = data.GetPtrFloat("jetRawEn");
    Int_t* jetID = NULL;
    //Int_t* jetID = data.GetPtrInt("jetID");

    //electron channel or muon channel
    isEle_ = 0;
    isMu_ = 0;
    z_mass_ = -99.;

    if ( outpath.Contains("JetHT")) {
      if ( (HLTJet >>10&1) != 1 && (HLTJet >>11&1) != 1 && (HLTJet >>12&1) != 1 && (HLTJet >>13&1) != 1
           && (HLTJet >>14&1) != 1 && (HLTJet >>15&1) != 1 && (HLTJet >>16&1) != 1 && (HLTJet >>17&1) != 1
           && (HLTJet >>18&1) != 1) continue;
    }

    /*
    TLorentzVector lep[10];
    //cout << "ele loop" << endl;
    //for electron channel
    vector<int> eleIndex;
    int ne = 0;
    for (int ie = 0; ie < nEle; ie++) {
      if (elePt[ie] < 15.) continue;
      if (fabs(eleSCEta[ie]) > 2.5) continue;
      if ( fabs(eleSCEta[ie]) > 1.4442 && fabs(eleSCEta[ie])<1.566 ) continue;
      //if ( Electron_SMZg(data,ie) == 0) continue;
      if ( ((eleID[ie] >> 3) &1) != 1 ) continue;

      lep[ne].SetPtEtaPhiM(elePt[ie], eleEta[ie], elePhi[ie], 0.5*0.001);

      //cout << "checking gen matching if mc" << endl;
      TLorentzVector mclep;
      if (data.HasMC()) {
	bool match = false;
	for (Int_t i=0; i<nMC; ++i) {
	  if (fabs(mcPID[i]) != 11) continue;
	  if ( ((mcStatusFlag[i] >> 0) & 1) == 0 ) continue;
	  mclep.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.5*0.001);
	  if ( ne == 0 && mclep.DeltaR(lep[0]) < 0.1 ) match = true;
	  if ( ne == 1 && mclep.DeltaR(lep[1]) < 0.1 ) match = true;
	}
	if ( !match) continue;
      }
      //cout << "push back ele index to vector " << endl;
      eleIndex.push_back(ie);
      ne++;
    }

    //cout << "after electron loop" << endl;
    if (ne >= 1) {
      //if (elePt[eleIndex[0]] < 25.) continue;
      isEle_ = 1;
      leptType_ = 11;
      TLorentzVector ll;
      ll = lep[0]+ lep[1];
      z_mass_ = ll.M();
    }

    //for muon channmul
    int nmu = 0;
    vector<int> muIndex;
    for (int imu = 0; imu < nMu; imu++) {
      if (muPt[imu] < 10.) continue;
      if (fabs(muEta[imu]) > 2.4) continue;
      //if ( Muon_SMZg(data,imu) == 0) continue;
      if ( ((muIDbit[imu] >> 2) &1) != 1) continue;
      if ( (muPFChIso[imu] + TMath::Max(0., muPFPhoIso[imu] + muPFNeuIso[imu] -0.5 * muPFPUIso[imu]))/muPt[imu] > 0.15) continue;

      lep[nmu].SetPtEtaPhiM(muPt[imu], muEta[imu], muPhi[imu], 0.5*0.001);

      TLorentzVector mclep;
      if (data.HasMC()) {
	bool match = false;
	for (Int_t i=0; i<nMC; ++i) {
	  if (fabs(mcPID[i]) != 13) continue;
	  if ( ((mcStatusFlag[i] >> 0) & 1) == 0 ) continue;
	  mclep.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.5*0.001);
	  if ( nmu == 0 && mclep.DeltaR(lep[0]) < 0.1 ) match = true;
	  if ( nmu == 1 && mclep.DeltaR(lep[1]) < 0.1 ) match = true;
	}
	if ( !match) continue;
      }
      muIndex.push_back(imu);
      nmu++;
    }

    //cout << "done muon loop" << endl;
    if (nmu >= 1) {
      //if (muPt[muIndex[0]] < 20.) continue;
      isMu_ = 1;
      leptType_ = 13;
      TLorentzVector ll;
      ll = lep[0]+ lep[1];
      z_mass_ = ll.M();
    }

    //comment out for background template
    if (ne <=1 && nmu <=1 ) continue;
    if (z_mass_ < 50.) continue;

    //for choosing FSR photon
    //if (z_mass_ < 40. || z_mass_ > 80.)  continue;
    //if (ne >=1 || nmu >=1 ) continue;
    */

    //if ( outpath.Contains("JetHT") )
    // if (nmu >=1 || ne >=1) continue; //for jet dataset                                                                                                        

    /************ select jet **************/
    TLorentzVector jet[10];
    nseljet = 0;
    if (outpath.Contains("JetHT") && nJet < 2) continue;

    //for JetHT dataset
    if ( outpath.Contains("JetHT") ){
      for (int ij = 0; ij < nJet; ij++) {
	if (jetPt[ij] < 45) continue;
	if (nseljet>5) continue;
	//if (jetPFLooseId[ij] != 1) continue;
	if ( jetID[ij] < 4) continue;
	jetPt_[nseljet] = jetPt[ij];
	jetPhi_[nseljet] = jetPhi[ij];
	jetEta_[nseljet] = jetEta[ij];

	jet[nseljet].SetPtEtaPhiE(jetPt[ij], jetEta[ij], jetPhi[ij], jetEn[ij]);
	nseljet++;
      }

      if ( nseljet < 2) continue;
    }

    //photon
    isEB_ = 0;
    isEE_ = 0;
    Int_t nSelPho = 0;

    genwei_ = 1;
    if (nPho < 1) continue;
    //cout << "photon loop" << endl;
    for (Int_t j=0; j<nPho; ++j) {
      if (PassPhotonPreselections(data, j) == 0) continue;
      pho.SetPtEtaPhiM(phoEt[j], phoEta[j], phoPhi[j], 0);

      if ( outpath.Contains("JetHT") ) {
        if (pho.DeltaR(jet[0]) < 0.7) continue;
        if (pho.DeltaR(jet[1]) < 0.7) continue;
      }

      //if (ne >=1 || nmu >=1)
      //if ( pho.DeltaR(lep[0]) < 0.7 || pho.DeltaR(lep[1]) < 0.7 ) continue; //remove FSR from lepton

      //for choosing FSR photon
      //if ( pho.DeltaR(lep[0]) > 0.8 && pho.DeltaR(lep[1]) > 0.8 ) continue; 
      //if ( pho.DeltaR(lep[0]) < 0.3 || pho.DeltaR(lep[1]) < 0.3 ) continue; 
      
      int genIndex = -1;
      if (data.HasMC()) {
	int genmatch = 0;
	int matchgenele = 0;
	int fragmentation = 0;
	
	for (Int_t i=0; i<nMC; ++i) {
	  bool match = false;
	  if (mcPID[i] != 22) continue;
	  if ( ((mcStatusFlag[i] >> 0) & 1) == 0 && ((mcStatusFlag[i] >> 1) & 1) == 0 ) continue;
	  if ( mcCalIsoDR03[i] < 5.) fragmentation = 1;
	  //if ( mcCalIsoDR03[i] > 5.) continue;
	  mcpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
	  if (pho.DeltaR(mcpho) < 0.2 && fabs(phoEt[j]-mcPt[i])/mcPt[i] < 0.2) {
	    //if (pho.DeltaR(mcpho) < 0.1) {
	    //genmatch = 1;
	    //genIndex = i;
	    //break;
	    match = true; 
	    genIndex = i;
	  }

	  if ( match && mcCalIsoDR03[genIndex] < 5.) {
	    genmatch = 1;
	    break;
	  }
	  else genmatch = 0;

	}
	
	//check matching gen ele
	for (Int_t i=0; i<nMC; ++i) {
          if (abs(mcPID[i]) != 11) continue;
          //if ( ((mcStatusFlag[i] >> 0) & 1) == 0 && ((mcStatusFlag[i] >> 1) & 1) == 0 ) continue;
          //if (((mcStatusFlag[i] >> 0) & 1) == 1 ) {
            mcpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
            if (pho.DeltaR(mcpho) < 0.2) {
              matchgenele = 1;
              break;
            }
	    //}
        }
	
	if (outpath.Contains("Zg_aMCatNLO") || outpath.Contains("WG")) {
	  if ( genmatch == 0) continue; //for signal
	  if ( mcCalIsoDR03[genIndex] > 5.) continue;
	}
	//if (outpath.Contains("DYJetsToLL") && genmatch == 1) continue; //for background
	if (outpath.Contains("DYJetsToLL") || outpath.Contains("QCD") || outpath.Contains("GJet")) {
	  if (genmatch == 1) continue;
	  //if ( fragmentation == 1) continue;
	  //if ( (fabs(phoSCEta[j]) < 1.5 && phoPFChIso[j] < 2.) || (fabs(phoSCEta[j]) > 1.5 && phoPFChIso[j] < 1.5) )
	  if (matchgenele == 1) continue;
        }
      }
      

      if (nSelPho == 0) {
	gamma_pt_ = phoEt[j];
	if (fabs(phoSCEta[j]) < 1.4442) isEB_ = 1;
	if (fabs(phoSCEta[j]) > 1.566 && fabs(phoSCEta[j]) < 2.5) isEE_ = 1;
	gamma_ssmva_ = PhotonSSMVA(data, j);
	gamma_ChIso_ = phoPFChIso[j];
	gamma_eta_ = phoEta[j];
	gamma_sceta_ = phoSCEta[j];
	gamma_phi_ = phoPhi[j];
	gamma_HoverE_ = phoHoverE[j];
	gamma_sigIetaIeta_ = phoSigmaIEtaIEtaFull5x5[j];
	gamma_ChWorstIso_ = phoPFChWorstIso[j];
	gamma_PhoIso_ = phoPFPhoIso[j];

	if (data.HasMC()) {
	  mcMomPID_ = mcMomPID[genIndex];
	  isogen_ = mcCalIsoDR03[genIndex];
	}

	//TLorentzVector mllg;
	//mllg = lep[0] + lep[1] + pho;
	//boss_mass_ = mllg.M();
	//boss_pt_ = mllg.Pt();

	outtree->Fill(); 
	nSelPho++;
      }
    }
    //cout<<nSelPho<<endl;
    
    counts++;
    
  } // event loop
  
  cout<<"processed "<<counts<<" events !"<<endl;
  fo->Write();
  fo->Close();
  delete fo;
}
