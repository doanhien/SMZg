#ifndef tree_h
#define tree_h

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <iostream>

float deltaPhi(float phi1, float phi2) {

  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return dPhi;

}

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return sqrt(pow(dEta,2) + pow (dPhi,2));

}

const int maxPar = 50;

Int_t run_;
Int_t lumi_;
Long64_t event_;
Int_t nVtx_;
Float_t vx_;
Float_t vy_;
Float_t vz_;
Float_t genWeight_;
Float_t puweigj;
Float_t puweigj_69p2nb;
Float_t puweigj_65nb;
Float_t puweigj_63nb;
Float_t weight;

Int_t aMCNLO;
Int_t nPU;
Int_t leptType;
Float_t lept0_pt;
Float_t lept0_eta;
Float_t lept0_phi;
Float_t lept0_sceta;

Float_t lept1_pt;
Float_t lept1_eta;
Float_t lept1_phi;
Float_t lept1_sceta;

Int_t   isEB;
Int_t   isEE;
Float_t gamma_SF;
Float_t gamma_pt;
Float_t gamma_eta;
Float_t gamma_sceta;
Float_t gamma_phi;
Float_t gamma_iso;
Float_t gamma_HoverE;
Float_t gamma_sigmaIetaIeta;
Float_t gamma_mva;
Float_t gamma_ssmva;
Float_t gamma_R9;
Float_t gamma_ChIso;
Float_t gamma_ChWorstIso;
Float_t gamma_EleVeto;

void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("weight",      &weight);
  tree->Branch("puweigj",     &puweigj);
  tree->Branch("puweigj_69p2nb",     &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",       &puweigj_65nb);
  tree->Branch("puweigj_63nb",       &puweigj_63nb);

  tree->Branch("leptType",                &leptType);
  tree->Branch("lept0_pt",                &lept0_pt);
  tree->Branch("lept0_eta",               &lept0_eta);
  tree->Branch("lept0_phi",               &lept0_phi);
  tree->Branch("lept0_sceta",             &lept0_sceta);

  tree->Branch("lept1_pt",                &lept1_pt);
  tree->Branch("lept1_eta",               &lept1_eta);
  tree->Branch("lept1_phi",               &lept1_phi);
  tree->Branch("lept1_sceta",             &lept1_sceta);

  tree->Branch("isEB",           &isEB);
  tree->Branch("isEE",           &isEE);
  tree->Branch("gamma_SF",       &gamma_SF);
  tree->Branch("gamma_pt",       &gamma_pt);
  tree->Branch("gamma_eta",      &gamma_eta);
  tree->Branch("gamma_sceta",    &gamma_sceta);
  tree->Branch("gamma_phi",      &gamma_phi);
  tree->Branch("gamma_iso",      &gamma_iso);
  tree->Branch("gamma_HoverE",   &gamma_HoverE);
  tree->Branch("gamma_sigmaIetaIeta", &gamma_sigmaIetaIeta);
  tree->Branch("gamma_mva",           &gamma_mva);
  tree->Branch("gamma_ssmva",         &gamma_ssmva);
  tree->Branch("gamma_R9",            &gamma_R9);
  tree->Branch("gamma_ChIso",         &gamma_ChIso);
  tree->Branch("gamma_ChWorstIso",    &gamma_ChWorstIso);
  tree->Branch("gamma_EleVeto",       &gamma_EleVeto);

}

Int_t     run;
Long64_t  event;
Int_t     lumis;
Bool_t    isData;
ULong64_t hlt;
ULong64_t hltPho;
Int_t nVtx;
Float_t vz;
Float_t vx; 
Float_t vy;
Float_t genWeight;
Float_t rho;

Int_t nEle ;
Int_t* eleCharge;
Int_t* eleChargeConsistent;
Float_t *eleEn;
Float_t* elePt ;
Float_t* eleEta ;
Float_t* elePhi ;
Float_t* eleHoverE ;
Float_t* eleEoverP ;
Float_t* eleEoverPInv ;
Float_t* eleSigmaIEtaIEta_Full5x5 ;
Float_t* eleSigmaIPhiIPhi;
Int_t* eleMissHits ;
Float_t* eleD0 ;
Float_t* eleDz ;
Float_t* eledEtaAtVtx ;
Float_t* eledPhiAtVtx ;
Float_t* eledEtaseedAtVtx ;
Int_t* eleConvVeto;
Int_t* eleEcalDrivenSeed;
Float_t* eleE1x5Full5x5;
Float_t* eleE2x5Full5x5;
Float_t* eleE5x5Full5x5;
Float_t* eleSCEta ;
Float_t* eleSCPhi;
Float_t* eleSCEtaWidth;
Float_t* eleSCPhiWidth;
Float_t* eleSCEn ;
Float_t* elePFChIso ;
Float_t* elePFPhoIso ;
Float_t* elePFNeuIso ;
Float_t* elePFPUIso ;
Float_t* elePFMiniIso ;
Short_t* eleID;
Float_t* eleIDMVA;
Float_t* eleIDMVAHZZ;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;
Float_t* eleGSFChi2;
Float_t* eleEcalEnErr;
Float_t* eleSIP;


Int_t    nMu;
Int_t*   muType;
Float_t* muPt;
Float_t* muEta;
Float_t* muPhi;
Int_t*   muCh;
Float_t*  muSIP;
Float_t* muChi2NDF;
Int_t*   muMuonHits;
Int_t*   muStations;
Int_t*   muTrkLayers;
Int_t*   muPixelHits;
Float_t* muInnerD0;
Float_t* muInnerDz;
Float_t* muD0;
Float_t* muDz;
Float_t* muBestTrkPtError;
Float_t* muBestTrkPt;
Int_t*   muBestTrkType;
Float_t* muPFChIso;
Float_t* muPFPhoIso;
Float_t* muPFNeuIso;
Float_t* muPFPUIso;
Float_t* muPFChIso03;
Float_t* muPFPhoIso03;
Float_t* muPFNeuIso03;
Float_t* muPFPUIso03;
Float_t* muPFMiniIso;
Float_t* muIsoTrk;
Short_t* muIDbit;


Int_t nPho;
Float_t* phoE;
Float_t* phoEt;
Float_t* phoSCEta;
Float_t* phoEta;
Float_t* phoPhi;
Float_t* phoCalibEt;
Float_t* phoHoverE;
Float_t* phoSigmaIEtaIEtaFull5x5;
Float_t* phoSigmaIEtaIPhiFull5x5;
Float_t* phoE2x2Full5x5;
Float_t* phoE5x5Full5x5;
Float_t* phoESEffSigmaRR;
Float_t* phoESEnP1;
Float_t* phoESEnP2;
Float_t* phoSCRawE;
Float_t* phoSCEtaWidth;
Float_t* phoSCPhiWidth;
Float_t* phoR9;
Float_t* phoPFChIso;
Float_t* phoPFNeuIso;
Float_t* phoPFPhoIso;
Float_t* phoPFChWorstIso;
Int_t* phohasPixelSeed;
Int_t* phoEleVeto;
Short_t* phoIDbit;
Float_t* phoIDMVA;
Float_t* phoEnErr;
Float_t* phoEn;


//gen level
Int_t nMC;
Int_t* mcPID;
Int_t* mcMomPID;
Int_t* mcGMomPID;
Float_t* mcPt;
Float_t* mcEt;
Float_t* mcPhi;
Float_t* mcEta;
Float_t* mcMomPt;
Float_t* mcMomEta;
UShort_t* mcStatusFlag;
Float_t* mcCalIsoDR03;


void readggtree(TreeReader &data) {

  run = data.GetInt("run");
  event = data.GetLong64("event");
  lumis = data.GetInt("lumis");
  hlt = data.GetLong64("HLTEleMuX");
  hltPho = data.GetLong64("HLTPho");
  nVtx = data.GetInt("nVtx");
  vx = data.GetFloat("vtx");
  vy = data.GetFloat("vty");
  vz = data.GetFloat("vtz");
  rho = data.GetFloat("rho");

  nEle = data.GetInt("nEle");
  eleCharge = data.GetPtrInt("eleCharge");
  elePt = data.GetPtrFloat("eleCalibPt");
  //elePt = data.GetPtrFloat("elePt");
  eleEn = data.GetPtrFloat("eleCalibEn");
  eleEta = data.GetPtrFloat("eleEta");
  elePhi = data.GetPtrFloat("elePhi");
  eleHoverE = data.GetPtrFloat("eleHoverE");
  eleEoverP = data.GetPtrFloat("eleEoverP");
  eleEoverPInv = data.GetPtrFloat("eleEoverPInv");
  eleSigmaIEtaIEta_Full5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
  eleMissHits = data.GetPtrInt("eleMissHits");
  eleD0 = data.GetPtrFloat("eleD0");
  eleDz = data.GetPtrFloat("eleDz");
  eledEtaAtVtx = data.GetPtrFloat("eledEtaAtVtx");
  eledPhiAtVtx = data.GetPtrFloat("eledPhiAtVtx");
  eledEtaseedAtVtx = data.GetPtrFloat("eledEtaseedAtVtx");
  eleConvVeto = data.GetPtrInt("eleConvVeto");
  eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  eleE1x5Full5x5 = data.GetPtrFloat("eleE1x5Full5x5");
  eleE2x5Full5x5 = data.GetPtrFloat("eleE2x5Full5x5");
  eleE5x5Full5x5 = data.GetPtrFloat("eleE5x5Full5x5");
  eleSCEta = data.GetPtrFloat("eleSCEta");
  eleSCPhi = data.GetPtrFloat("eleSCPhi");
  eleSCEtaWidth = data.GetPtrFloat("eleSCEtaWidth");
  eleSCPhiWidth = data.GetPtrFloat("eleSCPhiWidth");
  eleSCEn = data.GetPtrFloat("eleSCEn");
  elePFChIso = data.GetPtrFloat("elePFChIso");
  elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  elePFMiniIso = data.GetPtrFloat("elePFMiniIso");
  eleID = data.GetPtrShort("eleIDbit");
  eleIDMVA = data.GetPtrFloat("eleIDMVA");
  eleIDMVAHZZ = data.GetPtrFloat("eleIDMVAHZZ");
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");
  eleGSFChi2 = data.GetPtrFloat("eleGSFChi2");
  //eleEcalEnErr = data.GetPtrFloat("eleEcalEnErr");
  eleSIP = data.GetPtrFloat("eleSIP");

  nMu = data.GetInt("nMu");
  muType = data.GetPtrInt("muType");
  muPt = data.GetPtrFloat("muCalibPt");
  muEta = data.GetPtrFloat("muEta");
  muPhi = data.GetPtrFloat("muPhi");
  muCh = data.GetPtrInt("muCharge");
  muSIP = data.GetPtrFloat("muSIP");
  muChi2NDF = data.GetPtrFloat("muChi2NDF");
  muMuonHits = data.GetPtrInt("muMuonHits");
  muStations = data.GetPtrInt("muStations");
  muTrkLayers = data.GetPtrInt("muTrkLayers");
  muPixelHits = data.GetPtrInt("muPixelHits");
  muInnerD0 = data.GetPtrFloat("muInnerD0");
  muInnerDz = data.GetPtrFloat("muInnerDz");
  muD0 = data.GetPtrFloat("muD0");
  muDz = data.GetPtrFloat("muDz");
  muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError");
  muBestTrkPt = data.GetPtrFloat("muBestTrkPt");
  muBestTrkType = data.GetPtrInt("muBestTrkType");
  muPFChIso = data.GetPtrFloat("muPFChIso");
  muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
  muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
  muPFPUIso = data.GetPtrFloat("muPFPUIso");
  muPFChIso03 = data.GetPtrFloat("muPFChIso03");
  muPFNeuIso03 = data.GetPtrFloat("muPFNeuIso03");
  muPFPhoIso03 = data.GetPtrFloat("muPFPhoIso03");
  muPFPUIso03 = data.GetPtrFloat("muPFPUIso03");
  muPFMiniIso = data.GetPtrFloat("muPFMiniIso");
  muIsoTrk = data.GetPtrFloat("muIsoTrk");
  muIDbit = data.GetPtrShort("muIDbit");


  nPho = data.GetInt("nPho");
  phoE  = data.GetPtrFloat("phoE");
  phoEt = data.GetPtrFloat("phoEt");
  phoEta = data.GetPtrFloat("phoEta");
  phoSCEta = data.GetPtrFloat("phoSCEta");
  phoPhi = data.GetPtrFloat("phoPhi");
  phoCalibEt = data.GetPtrFloat("phoCalibEt");
  phoHoverE = data.GetPtrFloat("phoHoverE");
  phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
  phoE2x2Full5x5            = data.GetPtrFloat("phoE2x2Full5x5");
  phoE5x5Full5x5            = data.GetPtrFloat("phoE5x5Full5x5");
  phoESEnP1          = data.GetPtrFloat("phoESEnP1");
  phoESEnP2          = data.GetPtrFloat("phoESEnP2");
  phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");
  phoSCRawE         = data.GetPtrFloat("phoSCRawE");
  phoSCEtaWidth     = data.GetPtrFloat("phoSCEtaWidth");
  phoSCPhiWidth     = data.GetPtrFloat("phoSCPhiWidth");
  phoR9 = data.GetPtrFloat("phoR9");
  phoPFChIso = data.GetPtrFloat("phoPFChIso");
  phoPFNeuIso = data.GetPtrFloat("phoPFNeuIso");
  phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
  phoPFChWorstIso = data.GetPtrFloat("phoPFChWorstIso");
  phohasPixelSeed = data.GetPtrInt("phohasPixelSeed");
  phoEleVeto = data.GetPtrInt("phoEleVeto");
  phoIDbit = data.GetPtrShort("phoIDbit");
  phoIDMVA = data.GetPtrFloat("phoIDMVA");
  //phoEnErr = data.GetPtrFloat("phoEnErr");
  phoEn = data.GetPtrFloat("phoE");

  if (data.HasMC()) {
  nMC = data.GetInt("nMC");
  mcPID = data.GetPtrInt("mcPID");
  mcMomPID = data.GetPtrInt("mcMomPID");
  mcGMomPID = data.GetPtrInt("mcGMomPID");
  mcPt = data.GetPtrFloat("mcPt");
  mcEt = data.GetPtrFloat("mcEt");
  mcPhi = data.GetPtrFloat("mcPhi");
  mcEta = data.GetPtrFloat("mcEta");
  mcMomPt = data.GetPtrFloat("mcMomPt");
  mcMomEta = data.GetPtrFloat("mcMomEta");
  mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  genWeight = data.GetFloat("genWeight");
  mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");

}

}


#endif
