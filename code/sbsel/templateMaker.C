#include "TLorentzVector.h"
#include <vector>
#include <iostream>
#include <fstream>

#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TGraph.h"

#include "../ana/puweicalc.h"
#include "../ana/untuplizer.h"
#include "PhotonSelections.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"


using namespace std;

void templateMaker(const char** inpaths, int npaths, TString outpath = "minitree.root", Float_t xs=1, Float_t lumi=1, Int_t aMCatNLO=0) {

  /// inpaths  = array of paths to files with ggNtuplizer's TTrees;
  /// npaths   = size of inpaths;
  /// outpath  = path to output root file;
  /// xs       = cross section
  /// lumi     = luminosity
  /// aMCatNLO = 0 : not aMC@NLO, 1 : aMC@NLO (with negative weight) 

  // open tree(s) with events to be processed
  TreeReader data(inpaths, npaths);

  Long64_t event ;
  Float_t puwei_  = 1.; 
  Float_t mcwei_  = 1.;
  Float_t genwei_ = 1.;
  Float_t gamma_pt_, gamma_eta_, gamma_ssmva_, gamma_ChIso_;
  Int_t isEB_, isEE_;
  Int_t isEle_, isMu_;
  Int_t leptType_;
  Float_t z_mass_;
  Int_t   mcMomPID_; 
  Float_t isogen_;
  Float_t boss_mass_;
  Float_t boss_pt_;
  
  // prepare output tree
  TFile* fo = TFile::Open(outpath.Data(), "RECREATE");
  if (!fo || fo->IsZombie())  FATAL("TFile::Open() failed");
  //fo->mkdir("ggNtuplizer")->cd();
  fo->cd();
  TTree *outtreeMC_ = new TTree("outtree", "mini tree for building signal template");
  outtreeMC_->Branch("event",       &event,        "event/L");
  outtreeMC_->Branch("puwei",       &puwei_,       "puwei/F");
  outtreeMC_->Branch("mcwei",       &mcwei_,       "mcwei/F");
  outtreeMC_->Branch("genwei",      &genwei_,      "genwei/F");
  outtreeMC_->Branch("gamma_pt",    &gamma_pt_,    "gamma_pt/F"); 
  outtreeMC_->Branch("gamma_eta",   &gamma_eta_,   "gamma_eta/F"); 
  outtreeMC_->Branch("gamma_ssmva", &gamma_ssmva_, "gamma_ssmva/F");  
  outtreeMC_->Branch("gamma_ChIso", &gamma_ChIso_, "gamma_ChIso/F");
  outtreeMC_->Branch("isEB",        &isEB_,        "isEB/I"); 
  outtreeMC_->Branch("isEE",        &isEE_,        "isEE/I"); 
  outtreeMC_->Branch("isEle",       &isEle_,       "isEle/I"); 
  outtreeMC_->Branch("isMu",        &isMu_,        "isMu/I"); 
  outtreeMC_->Branch("leptType",    &leptType_,    "leptType/I");
  outtreeMC_->Branch("z_mass",      &z_mass_,      "z_mass/F");
  outtreeMC_->Branch("mcMomPID",    &mcMomPID_,    "mcMomPID/I");
  outtreeMC_->Branch("isogen",      &isogen_,      "isogen/F"); 
  outtreeMC_->Branch("boss_mass",   &boss_mass_,   "boss_mass/F");
  outtreeMC_->Branch("boss_pt",     &boss_pt_,     "boss_pt/F");



  // pileup reweighting for MC
  PUWeightCalculator puCalc;
  if (data.HasMC()) puCalc.Init("../ana/external/PU_histo_13TeV_GoldenJSON_65550nb.root");

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
  Int_t nsellep = 0;

  //ofstream fout;
  //fout.open("event_mcPID.txt");

  //TFile* fss_rw = TFile::Open("external/transformation_Moriond17_AfterPreApr_v1.root");
  TFile* fss_rw = TFile::Open("../showershapeCorr/transformation_pho_presel_BDTUpto6000.root");
  TGraph *tgr[20];
  tgr[0] = (TGraph*) fss_rw->Get("transfEtaWidthEB");
  tgr[1] = (TGraph*) fss_rw->Get("transfS4EB");
  tgr[2] = (TGraph*) fss_rw->Get("transffull5x5R9EB");
  tgr[3] = (TGraph*) fss_rw->Get("transffull5x5sieieEB");
  tgr[4] = (TGraph*) fss_rw->Get("transffull5x5sieipEB");
  tgr[5] = (TGraph*) fss_rw->Get("transfrhoEB");
  tgr[6] = (TGraph*) fss_rw->Get("transfPhiWidthEB");


  tgr[7] = (TGraph*) fss_rw->Get("transfEtaWidthEE");
  tgr[8] = (TGraph*) fss_rw->Get("transfS4EE");
  tgr[9] = (TGraph*) fss_rw->Get("transffull5x5R9EE");
  tgr[10] = (TGraph*) fss_rw->Get("transffull5x5sieieEE");
  tgr[11] = (TGraph*) fss_rw->Get("transffull5x5sieipEE");
  tgr[12] = (TGraph*) fss_rw->Get("transfrhoEE");
  tgr[13] = (TGraph*) fss_rw->Get("transfPhiWidthEE");


  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    //for (Long64_t ev = 0; ev < 50000000; ev++) {
    if ((ev % 50000 ) == 0) cout << "processing event " << ev+1 << " in " << data.GetEntriesFast() << endl;
    data.GetEntry(ev);
    
    Int_t run_  = data.GetInt("run");
    Float_t rho = data.GetFloat("rho");
    Int_t lumis_ = data.GetInt("lumis");
    Long64_t event_ = data.GetLong64("event");
    //nVtx_  = data.GetInt("nVtx");
    
    event = event_;
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
    float* mcPt     = NULL;
    float* mcEta    = NULL;
    float* mcPhi    = NULL;
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
    float* phoEt      = data.GetPtrFloat("phoEt");
    float* phoSCEta   = data.GetPtrFloat("phoSCEta");
    float* phoEta     = data.GetPtrFloat("phoEta");
    float* phoPhi     = data.GetPtrFloat("phoPhi");
    float* phoPFChIso = data.GetPtrFloat("phoPFChIso");   

    Int_t nEle        = data.GetInt("nEle");
    Float_t* elePt    = data.GetPtrFloat("elePt");
    Float_t* eleEta   = data.GetPtrFloat("eleEta");
    Float_t* eleSCEta = data.GetPtrFloat("eleSCEta");
    Float_t* elePhi   = data.GetPtrFloat("elePhi");

    Int_t nMu         = data.GetInt("nMu");
    Float_t* muPt     = data.GetPtrFloat("muPt");
    Float_t* muEta    = data.GetPtrFloat("muEta");
    Float_t* muPhi    = data.GetPtrFloat("muPhi");

    //electron channel or muon channel
    isEle_ = 0;
    isMu_ = 0;
    z_mass_ = -99.;


    TLorentzVector lep[10];
    //cout << "ele loop" << endl;
    //for electron channel
    vector<int> eleIndex;
    vector<int> selectedEle;
    //Electron_SMZg(data, selectedEle);
    passEleId(data, selectedEle);

    int ne = 0;
    if ( selectedEle.size() >= 1 ) {
      for (vector<int>::const_iterator ie = selectedEle.begin(); ie != selectedEle.end(); ie++) {
	if (data.HasMC()) {
	  if ( eleMatcher(data, *ie)!=1 ) continue;
	}
	if (elePt[*ie] < 15.) continue;

	eleIndex.push_back(*ie);
	lep[ne].SetPtEtaPhiM(elePt[*ie],eleEta[*ie],elePhi[*ie], 0.511*0.001);

	ne++;
      }
    } //done electron loop


    if (ne >= 1) {
      if (elePt[eleIndex[0]] < 25.) continue;
      isEle_ = 1;
      leptType_ = 11;
      TLorentzVector ll;
      ll = lep[0]+ lep[1];
      z_mass_ = ll.M();
    }

    //for muon channmul
    int nmu = 0;
    vector<int> muIndex;
    vector<int> selectedMu;
    //Muon_SMZg(data, selectedMu);
    passMuonId(data, selectedMu);

    if ( selectedMu.size() >= 1 ) {
      for (vector<int>::const_iterator imu = selectedMu.begin(); imu != selectedMu.end(); imu++) {
	if (data.HasMC() ) {
	  if (muonMatcher(data, *imu) !=1 ) continue;
	}
	muIndex.push_back(*imu);
	lep[nmu].SetPtEtaPhiM(muPt[*imu], muEta[*imu], muPhi[*imu], 0.104);
	nmu++;
      }
    }

    //cout << "done muon loop" << endl;
    if (nmu >= 1) {
      if (muPt[muIndex[0]] < 20.) continue;
      isMu_ = 1;
      leptType_ = 13;
      TLorentzVector ll;
      ll = lep[0]+ lep[1];
      z_mass_ = ll.M();
    }

    //comment out for background template
    //if (ne < 1 && nmu < 1 ) continue;
    if (ne >= 1 && nmu >= 1 ) nsellep++;
    //if (z_mass_ < 50.) continue;

    //for choosing FSR photon
    //if (z_mass_ < 40. || z_mass_ > 80.)  continue;
    //if (ne >=1 || nmu >=1 ) continue;


    //photon
    isEB_ = 0;
    isEE_ = 0;
    Int_t nSelPho = 0;
    /*
    ifstream infile;
    infile.open("event.txt");
    Long64_t Ev;
    ///int id, momid, status0, status1;
    bool chosedEv = false;
    if (infile.is_open()) {
      while (infile >> Ev) {
	if (Ev != event_) continue;
	else chosedEv = true;
      }
    }
    infile.close();
    */
    if (nPho < 1) continue;
    genwei_ = 1;
    for (Int_t j=0; j<nPho; ++j) {
      if (phoEt[j] < 15.) continue;
      if (PassPhotonPreselections(data, j) == 0) continue;

      pho.SetPtEtaPhiM(phoEt[j], phoEta[j], phoPhi[j], 0);
      

      if ( (ne == 1 || nmu == 1) &&  pho.DeltaR(lep[0]) < 0.7 ) continue;
      else if (ne >1 || nmu >1)
	if ( pho.DeltaR(lep[0]) < 0.7 || pho.DeltaR(lep[1]) < 0.7 ) continue; //remove FSR from lepton

      //for choosing FSR photon
      //if ( pho.DeltaR(lep[0]) > 0.8 && pho.DeltaR(lep[1]) > 0.8 ) continue; 
      //if ( pho.DeltaR(lep[0]) < 0.3 || pho.DeltaR(lep[1]) < 0.3 ) continue; 
      

      int genIndex = -1;
      if (data.HasMC()) {
	if ( outpath.Contains("Zg") || outpath.Contains("Wg") || outpath.Contains("GJet")) { 
	  if ( phoMatcher_prompt(data, j) != 1) continue;
	} 
	else {      
	  if ( phoMatcher_prompt(data, j) ) continue;
	  //if (pho_matchele(data, j) ) continue;
	  //if (pho_matchfsr(data, j) ) continue;

	}
      }

      
      if (nSelPho == 0) {
	gamma_pt_ = phoEt[j];
	gamma_eta_ = phoEta[j];
	if (fabs(phoSCEta[j]) < 1.4442) isEB_ = 1;
	if (fabs(phoSCEta[j]) > 1.566 && fabs(phoSCEta[j]) < 2.5) isEE_ = 1;
	gamma_ssmva_ = PhotonSSMVA(data, j, tgr);

	//rho-correction charged isolation
	float chEA[] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
	int bin;
	if (phoSCEta[j] < 1.0) bin = 0;
	else if (phoSCEta[j] >= 1.0 && phoSCEta[j] < 1.479) bin = 1;
	else if (phoSCEta[j] >= 1.479 && phoSCEta[j] < 2.0) bin = 2;
	else if (phoSCEta[j] >= 2.0 && phoSCEta[j] < 2.2) bin = 3;
	else if (phoSCEta[j] >= 2.2 && phoSCEta[j] < 2.3) bin = 4;
	else if (phoSCEta[j] >= 2.3 && phoSCEta[j] < 2.4) bin = 5;
	else bin = 6; // phoSCEta[j] >= 2.4

	double corr = phoPFChIso[j] - chEA[bin] * rho;
	double corrchIso = TMath::Max(corr, 0.);

	//gamma_ChIso_ = phoPFChIso[j];
	gamma_ChIso_ = corrchIso;
	if (data.HasMC()) {
	  mcMomPID_ = mcMomPID[genIndex];
	  isogen_ = mcCalIsoDR03[genIndex];
	  //if ( gamma_ssmva_>0.7 && gamma_ssmva_<0.9 && phoEt[j]>15 && phoEt[j]<25 && fabs(phoSCEta[j]) < 1.4442 && gamma_ChIso_<2.)
	    //fout << event_ << "\t" << gamma_pt_ << "\t" << phoSCEta[j] << "\t" << phoEta[j] << endl;

	  //if ( chosedEv)
	  //cout << "event: " << event_ << "\t match_prompt: " << phoMatcher_prompt(data, j) 
	  //	 << "\t match_ele: " << pho_matchele(data, j) 
	  //	 << "\t match fsr: " << pho_matchfsr(data, j) << endl;

	}

	TLorentzVector mllg;
	//mllg = lep[0] + lep[1] + pho;
	boss_mass_ = mllg.M();
	boss_pt_ = mllg.Pt();

	outtreeMC_->Fill(); 
	nSelPho++;
      }
    }
    //infile.close(); 

   
    counts++;
    
  } // event loop
  
  cout<<"processed "<<counts<<" events !"<<endl;
  cout << "number of event having >= 1 ele or muon: " << nsellep << endl;
  //fout.close();
  fo->Write();
  fo->Close();
  delete fo;
}
