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
Float_t puweigj_69nb;
Float_t puweigj_69p2nb;
Float_t puweigj_65nb;
Float_t puweigj_63nb;

Int_t     trig_doublePho70;
Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele35_WPTight;
Int_t     trig_Ele17_Ele12;
Int_t     trig_Ele23_Ele12;
Int_t     trig_DoubleEle33;

Int_t     ne;
Int_t     nTag;
Int_t     TagIndex;
Float_t   tag_elePt;
Float_t   tag_eleEta;
Float_t   tag_elePhi;
Float_t   tag_eleSCEta;
Int_t     tag_ele_MatchTrig27;
Int_t     tag_ele_MatchTrig35;
Int_t     nProbe;
Int_t     ProbeIndex;
Float_t   probe_elePt;
Float_t   probe_eleEta;
Float_t   probe_elePhi;
Float_t   probe_eleSCEta;
Int_t     probe_ele_MatchTrig_Leg1_23;
Float_t   Zm;
Float_t   Zeta;
Float_t   Zy;
Int_t     npair;

Int_t     nProbe_matchPho;
Float_t   probe_phoPt;
Float_t   probe_phoEta;
Float_t   probe_phoPhi;
Float_t   probe_phoSCEta;
Float_t   probe_phoChIso;
Int_t     probe_phoEleVeto;
Float_t   Zmass_;
Float_t   Zeta_;
Float_t   Zy_;
Int_t     npaireg;


void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj_69nb",     &puweigj_69nb);
  tree->Branch("puweigj_69p2nb",   &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",     &puweigj_65nb);
  tree->Branch("puweigj_63nb",     &puweigj_63nb);

  tree->Branch("trig_doublePho70",        &trig_doublePho70);
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele35_WPTight",      &trig_Ele35_WPTight);
  tree->Branch("trig_Ele17_Ele12",        &trig_Ele17_Ele12);
  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);

  tree->Branch("ne",           &ne);
  tree->Branch("nTag",         &nTag);
  tree->Branch("TagIndex",     &TagIndex);
  tree->Branch("tag_elePt",       &tag_elePt);
  tree->Branch("tag_eleEta",      &tag_eleEta);
  tree->Branch("tag_elePhi",      &tag_elePhi);
  tree->Branch("tag_eleSCEta",    &tag_eleSCEta);
  tree->Branch("tag_ele_MatchTrig27", &tag_ele_MatchTrig27);
  tree->Branch("tag_ele_MatchTrig35", &tag_ele_MatchTrig35);
  tree->Branch("nProbe",       &nProbe);
  tree->Branch("ProbeIndex",   &ProbeIndex);
  tree->Branch("probe_elePt",       &probe_elePt);
  tree->Branch("probe_eleEta",      &probe_eleEta);
  tree->Branch("probe_elePhi",      &probe_elePhi);
  tree->Branch("probe_eleSCEta",    &probe_eleSCEta);
  tree->Branch("probe_ele_MatchTrig_Leg1_23",  &probe_ele_MatchTrig_Leg1_23);
  tree->Branch("Zm",           &Zm);
  tree->Branch("Zeta",         &Zeta);
  tree->Branch("Zy",           &Zy);
  tree->Branch("npair",        &npair);

  tree->Branch("nProbe_matchPho",   &nProbe_matchPho);
  tree->Branch("probe_phoPt",       &probe_phoPt);
  tree->Branch("probe_phoEta",      &probe_phoEta);
  tree->Branch("probe_phoPhi",      &probe_phoPhi);
  tree->Branch("probe_phoSCEta",    &probe_phoSCEta);
  tree->Branch("probe_phoEleVeto",  &probe_phoEleVeto);
  tree->Branch("probe_phoChIso",    &probe_phoChIso);

  tree->Branch("Zmass",        &Zmass_);
  tree->Branch("Zeta",         &Zeta_);
  tree->Branch("Zy",           &Zy_);
  tree->Branch("npaireg",      &npaireg);

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
Float_t* eleIDMVANonTrg;
Int_t*   eleFiredSingleTrgs;
//Long64_t*   eleFiredSingleTrgs;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;


Int_t    nMu;
Float_t* muPt;
Float_t* muEta;
Float_t* muPhi;
Int_t*   muCh;
Float_t* muChi2NDF;
Int_t*   muMuonHits;
Int_t*   muStations;
Int_t*   muTrkLayers;
Int_t*   muPixelHits;
Float_t* muInnerD0;
Float_t* muInnerDz;
Float_t* muPFChIso;
Float_t* muPFPhoIso;
Float_t* muPFNeuIso;
Float_t* muPFPUIso;

Int_t nPho;
Float_t* phoEt;
Float_t* phoSCEta;
Float_t* phoSCPhi;
Float_t* phoEta;
Float_t* phoPhi;
Float_t* phoCalibEt;
Float_t* phoHoverE;
Float_t* phoSigmaIEtaIEtaFull5x5;
Float_t* phoR9;
Float_t* phoPFChIso;
Float_t* phoPFNeuIso;
Float_t* phoPFPhoIso;
Int_t* phohasPixelSeed;
Int_t* phoEleVeto;
Short_t* phoIDbit;
Float_t* phoIDMVA;
Float_t* phoPFChWorstIso;
Int_t* phoFiredDoubleTrgs;
//Long64_t* phoFiredDoubleTrgs;



//gen level
Int_t nMC;
Int_t* mcPID;
Int_t* mcMomPID;
Float_t* mcPt;
Float_t* mcPhi;
Float_t* mcEta;
UShort_t* mcStatusFlag;


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
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  eleFiredSingleTrgs = data.GetPtrInt("eleFiredSingleTrgs");
  //eleFiredSingleTrgs = data.GetPtrLong64("eleFiredSingleTrgs");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");

  nMu = data.GetInt("nMu");
  muPt = data.GetPtrFloat("muPt");
  muEta = data.GetPtrFloat("muEta");
  muPhi = data.GetPtrFloat("muPhi");
  muCh = data.GetPtrInt("muCharge");
  muChi2NDF = data.GetPtrFloat("muChi2NDF");
  muMuonHits = data.GetPtrInt("muMuonHits");
  muStations = data.GetPtrInt("muStations");
  muTrkLayers = data.GetPtrInt("muTrkLayers");
  muPixelHits = data.GetPtrInt("muPixelHits");
  muInnerD0 = data.GetPtrFloat("muInnerD0");
  muInnerDz = data.GetPtrFloat("muInnerDz");
  muPFChIso = data.GetPtrFloat("muPFChIso");
  muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
  muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
  muPFPUIso = data.GetPtrFloat("muPFPUIso");

  nPho = data.GetInt("nPho");
  phoEt = data.GetPtrFloat("phoEt");
  phoEta = data.GetPtrFloat("phoEta");
  phoPhi = data.GetPtrFloat("phoPhi");
  phoSCEta = data.GetPtrFloat("phoSCEta");
  phoSCPhi = data.GetPtrFloat("phoSCPhi");
  phoCalibEt = data.GetPtrFloat("phoCalibEt");
  //phoCalibEt = data.GetPtrFloat("phoEt");
  phoHoverE = data.GetPtrFloat("phoHoverE");
  phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  phoR9 = data.GetPtrFloat("phoR9");
  phoPFChIso = data.GetPtrFloat("phoPFChIso");
  phoPFNeuIso = data.GetPtrFloat("phoPFNeuIso");
  phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
  phoPFChWorstIso = data.GetPtrFloat("phoPFChWorstIso");
  phohasPixelSeed = data.GetPtrInt("phohasPixelSeed");
  phoEleVeto = data.GetPtrInt("phoEleVeto");
  phoIDbit = data.GetPtrShort("phoIDbit");
  phoIDMVA = data.GetPtrFloat("phoIDMVA");
  phoFiredDoubleTrgs = data.GetPtrInt("phoFiredDoubleTrgs");
  //phoFiredDoubleTrgs = data.GetPtrLong64("phoFiredDoubleTrgs");

  if (data.HasMC()) {
  nMC = data.GetInt("nMC");
  mcPID = data.GetPtrInt("mcPID");
  mcMomPID = data.GetPtrInt("mcMomPID");
  mcPt = data.GetPtrFloat("mcPt");
  mcPhi = data.GetPtrFloat("mcPhi");
  mcEta = data.GetPtrFloat("mcEta");
  mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  genWeight = data.GetFloat("genWeight");
}

}


#endif
