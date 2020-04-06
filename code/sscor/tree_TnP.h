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
Int_t   isPVGood_;
Float_t rho_;

Int_t     trig_doublePho60;
Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele27_Eta2p1_WPTight;
Int_t     trig_Ele17_Ele12;
Int_t     trig_Ele23_Ele12;
Int_t     trig_DoubleEle33;

Float_t   pfMet;
Float_t   pfMetPhi;
Int_t     ne;
Int_t     nTag;
Int_t     TagIndex;
Float_t   Tag_Pt;
Float_t   Tag_Eta;
Float_t   Tag_Phi;
Float_t   Tag_SCEta;
Float_t   Tag_sieie;
Float_t   Tag_chIso;
Float_t   Tag_neuIso;
Float_t   Tag_phoIso;
Float_t   Tag_hoe;
Float_t   Tag_eop;
Float_t   Tag_ch;
Float_t   Tag_dEtaSeed;
Float_t   Tag_dEtaVtx;
Float_t   Tag_dPhiVtx;
Float_t   Tag_IDMVA;
Int_t     ele1_MatchTrig27;
Float_t   mW_Tag_Met;
Float_t   Tag_RecoSF;
Float_t   Tag_SelSF;

Int_t     nProbe;
Int_t     ProbeIndex;
Int_t     nPassId;
Int_t     nFailId;
Int_t     nPassIdIso;
Int_t     nFailIdIso;
Float_t   Probe_Pt;
Float_t   Probe_Eta;
Float_t   Probe_Phi;
Float_t   Probe_SCEta;
Float_t   Probe_sieie;
Float_t   Probe_chIso;
Float_t   Probe_neuIso;
Float_t   Probe_phoIso;
Float_t   Probe_hoe;
Float_t   Probe_eop;
Float_t   Probe_ch;
Float_t   Probe_dEtaSeed;
Float_t   Probe_dEtaVtx;
Float_t   Probe_dPhiVtx;
Float_t   Probe_miniIso;
Int_t     LooseID;
Int_t     MediumID;
Int_t     TightID;
Int_t     LooseIDOnly;
Float_t   Probe_PhoEleVeto;
Float_t   Probe_PhoChIso;
Int_t     passPhoID_Zg;
Int_t     pho_presel;
Float_t   Probe_ssmva;
Float_t   Probe_phoSF;
Float_t   Probe_ChWorstIso;
Float_t   Probe_SCEtaWidth;
Float_t   Probe_s4Full5x5 ;
Float_t   Probe_R9;
Float_t   Probe_sieip;
Float_t   Probe_SCRawE;
Float_t   Probe_scphiwidth;
Float_t   Probe_ssmva_rw;
Float_t   Probe_SCEtaWidth_rw;
Float_t   Probe_s4Full5x5_rw ;
Float_t   Probe_R9_rw        ;
Float_t   Probe_sieie_rw;
Float_t   Probe_sieip_rw;
Float_t   rho_rw;
Float_t   Probe_SCPhiWidth_rw;
Float_t   diff_mva;

Float_t   Zm;
Float_t   Zeta;
Float_t   Zpt;

void inittree(TTree* tree) {

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("nVtx",        &nVtx_);
  tree->Branch("vx",          &vx_);
  tree->Branch("vy",          &vy_);
  tree->Branch("vz",          &vz_);
  tree->Branch("isPVGood",    &isPVGood_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj_69nb",     &puweigj_69nb);
  tree->Branch("puweigj_69p2nb",   &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",     &puweigj_65nb);
  tree->Branch("puweigj_63nb",     &puweigj_63nb);
  tree->Branch("rho",              &rho_);

  tree->Branch("trig_doublePho60",        &trig_doublePho60);
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele27_Eta2p1_WPTight",      &trig_Ele27_Eta2p1_WPTight);
  tree->Branch("trig_Ele17_Ele12",        &trig_Ele17_Ele12);
  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);

  tree->Branch("pfMet",        &pfMet);
  tree->Branch("pfMetPhi",     &pfMetPhi);
  tree->Branch("ne",           &ne);
  tree->Branch("nTag",         &nTag);
  tree->Branch("TagIndex",     &TagIndex);
  tree->Branch("Tag_Pt",       &Tag_Pt);
  tree->Branch("Tag_Eta",      &Tag_Eta);
  tree->Branch("Tag_Phi",      &Tag_Phi);
  tree->Branch("Tag_SCEta",    &Tag_SCEta);
  tree->Branch("Tag_sieie",    &Tag_sieie);
  tree->Branch("Tag_chIso",    &Tag_chIso);
  tree->Branch("Tag_neuIso",   &Tag_neuIso);
  tree->Branch("Tag_phoIso",   &Tag_phoIso);
  tree->Branch("Tag_hoe",      &Tag_hoe);
  tree->Branch("Tag_eop",      &Tag_eop);
  tree->Branch("Tag_ch",       &Tag_ch);
  tree->Branch("Tag_dEtaSeed",     &Tag_dEtaSeed);
  tree->Branch("Tag_dEtaVtx",      &Tag_dEtaVtx);
  tree->Branch("Tag_dPhiVtx",      &Tag_dPhiVtx);
  tree->Branch("Tag_IDMVA",        &Tag_IDMVA);
  tree->Branch("ele1_MatchTrig27", &ele1_MatchTrig27);
  tree->Branch("mW_Tag_Met",       &mW_Tag_Met);
  tree->Branch("Tag_RecoSF",       &Tag_RecoSF);
  tree->Branch("Tag_SelSF",        &Tag_SelSF);

  tree->Branch("nProbe",         &nProbe);
  tree->Branch("nPassId",        &nPassId);
  tree->Branch("nFailId",        &nFailId);
  tree->Branch("nPassIdIso",     &nPassIdIso);
  tree->Branch("nFailIdIso",     &nFailIdIso);
  tree->Branch("ProbeIndex",     &ProbeIndex);
  tree->Branch("Probe_Pt",       &Probe_Pt);
  tree->Branch("Probe_Eta",      &Probe_Eta);
  tree->Branch("Probe_Phi",      &Probe_Phi);
  tree->Branch("Probe_SCEta",    &Probe_SCEta);
  tree->Branch("Probe_chIso",    &Probe_chIso);
  tree->Branch("Probe_neuIso",   &Probe_neuIso);
  tree->Branch("Probe_phoIso",   &Probe_phoIso);
  tree->Branch("Probe_hoe",      &Probe_hoe);
  tree->Branch("Probe_eop",      &Probe_eop);
  tree->Branch("Probe_dEtaSeed", &Probe_dEtaSeed);
  tree->Branch("Probe_dEtaVtx",  &Probe_dEtaVtx);
  tree->Branch("Probe_dPhiVtx",  &Probe_dPhiVtx);
  tree->Branch("Probe_ch",       &Probe_ch);
  tree->Branch("Probe_miniIso",  &Probe_miniIso);
  tree->Branch("LooseID",        &LooseID);
  tree->Branch("MediumID",       &MediumID);
  tree->Branch("TightID",        &TightID);
  tree->Branch("LooseIDOnly",    &LooseIDOnly);
  tree->Branch("Probe_PhoEleVeto",  &Probe_PhoEleVeto);
  tree->Branch("Probe_PhoChIso",    &Probe_PhoChIso);
  tree->Branch("passPhoID_Zg",      &passPhoID_Zg);
  tree->Branch("pho_presel",        &pho_presel);
  tree->Branch("Probe_ssmva",       &Probe_ssmva);
  tree->Branch("Probe_phoSF",       &Probe_phoSF);
  tree->Branch("Probe_ChWorstIso",  &Probe_ChWorstIso);
  tree->Branch("Probe_sieip",       &Probe_sieip);
  tree->Branch("Probe_SCRawE",      &Probe_SCRawE);
  tree->Branch("Probe_scphiwidth",  &Probe_scphiwidth);

  tree->Branch("Probe_SCEtaWidth",  &Probe_SCEtaWidth);
  tree->Branch("Probe_s4Full5x5",   &Probe_s4Full5x5);
  tree->Branch("Probe_R9",          &Probe_R9);
  tree->Branch("Probe_sieie",       &Probe_sieie);

  tree->Branch("Probe_ssmva_rw",    &Probe_ssmva_rw);
  tree->Branch("diff_mva",          &diff_mva);
  tree->Branch("Probe_SCEtaWidth_rw", &Probe_SCEtaWidth_rw);
  tree->Branch("Probe_s4Full5x5_rw",  &Probe_s4Full5x5_rw);
  tree->Branch("Probe_R9_rw",         &Probe_R9_rw);
  tree->Branch("Probe_sieie_rw",      &Probe_sieie_rw);
  tree->Branch("Probe_sieip_rw",      &Probe_sieip_rw );
  tree->Branch("rho_rw",              &rho_rw);
  tree->Branch("Probe_SCPhiWidth_rw", &Probe_SCPhiWidth_rw);

  tree->Branch("Zm",             &Zm);
  tree->Branch("Zeta",           &Zeta);
  tree->Branch("Zpt",            &Zpt);

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
Int_t   isPVGood;
Float_t genWeight;
Float_t rho;
Float_t pfMET;
Float_t pfMETPhi;

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
Long64_t*   eleFiredSingleTrgs;
//Int_t*   eleFiredSingleTrgs;
Float_t* elePFClusEcalIso;
Float_t* elePFClusHcalIso;
Float_t* eleDr03TkSumPt;

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
Float_t* phoPFChWorstIso;
Int_t* phohasPixelSeed;
Int_t* phoEleVeto;
Short_t* phoIDbit;
Float_t* phoIDMVA;
Float_t* phoSCEtaWidth;
Float_t* phoSCPhiWidth;
Float_t* phoE2x2Full5x5;
Float_t* phoE5x5Full5x5;
Float_t* phoSigmaIEtaIPhiFull5x5;
Float_t* phoSCRawE;




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
  isPVGood = data.GetBool("isPVGood");
  rho = data.GetFloat("rho");
  pfMET = data.GetFloat("pfMET");
  pfMETPhi = data.GetFloat("pfMETPhi");
  nEle = data.GetInt("nEle");
  eleCharge = data.GetPtrInt("eleCharge");
  //elePt = data.GetPtrFloat("eleCalibPt");
  elePt = data.GetPtrFloat("elePt");
  //eleEn = data.GetPtrFloat("eleCalibEn");
  eleEn = data.GetPtrFloat("eleEn");
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
  //eledEtaseedAtVtx = data.GetPtrFloat("eledEtaseedAtVtx");
  //eleConvVeto = data.GetPtrInt("eleConvVeto");
  //eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
  //eleE1x5Full5x5 = data.GetPtrFloat("eleE1x5Full5x5");
  //eleE2x5Full5x5 = data.GetPtrFloat("eleE2x5Full5x5");
  //eleE5x5Full5x5 = data.GetPtrFloat("eleE5x5Full5x5");
  eleSCEta = data.GetPtrFloat("eleSCEta");
  eleSCPhi = data.GetPtrFloat("eleSCPhi");
  eleSCEtaWidth = data.GetPtrFloat("eleSCEtaWidth");
  eleSCPhiWidth = data.GetPtrFloat("eleSCPhiWidth");
  eleSCEn = data.GetPtrFloat("eleSCEn");
  elePFChIso = data.GetPtrFloat("elePFChIso");
  elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  //elePFMiniIso = data.GetPtrFloat("elePFMiniIso");
  eleID = data.GetPtrShort("eleIDbit");
  //eleIDMVA = data.GetPtrFloat("eleIDMVA");
  //eleIDMVANonTrg = data.GetPtrFloat("eleIDMVANonTrg");
  //eleFiredSingleTrgs = data.GetPtrInt("eleFiredSingleTrgs");
  eleFiredSingleTrgs = data.GetPtrLong64("eleFiredSingleTrgs");
  elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
  elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
  //eleDr03TkSumPt = data.GetPtrFloat("eleDr03TkSumPt");

  nPho = data.GetInt("nPho");
  //phoEt = data.GetPtrFloat("phoEt");
  phoEt = data.GetPtrFloat("phoCalibEt");
  phoEta = data.GetPtrFloat("phoEta");
  phoPhi = data.GetPtrFloat("phoPhi");
  phoSCEta = data.GetPtrFloat("phoSCEta");
  phoSCPhi = data.GetPtrFloat("phoSCPhi");
  phoCalibEt = data.GetPtrFloat("phoCalibEt");
  phoHoverE = data.GetPtrFloat("phoHoverE");
  phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  phoR9 = data.GetPtrFloat("phoR9");
  phoSCRawE         = data.GetPtrFloat("phoSCRawE");
  phoSCEtaWidth     = data.GetPtrFloat("phoSCEtaWidth");
  phoSCPhiWidth     = data.GetPtrFloat("phoSCPhiWidth");
  phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
  phoE2x2Full5x5            = data.GetPtrFloat("phoE2x2Full5x5");
  phoE5x5Full5x5            = data.GetPtrFloat("phoE5x5Full5x5");
  phoPFChIso = data.GetPtrFloat("phoPFChIso");
  phoPFNeuIso = data.GetPtrFloat("phoPFNeuIso");
  phoPFPhoIso = data.GetPtrFloat("phoPFPhoIso");
  phoPFChWorstIso = data.GetPtrFloat("phoPFChWorstIso");
  phohasPixelSeed = data.GetPtrInt("phohasPixelSeed");
  phoEleVeto = data.GetPtrInt("phoEleVeto");
  phoIDbit = data.GetPtrShort("phoIDbit");
  phoIDMVA = data.GetPtrFloat("phoIDMVA");

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
