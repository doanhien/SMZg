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
Float_t rho_;

Int_t     trig_doublePho60;
Int_t     trig_TriPho303010;
Int_t     trig_TriPho35355;
Int_t     trig_Ele27_WPTight;
Int_t     trig_Ele17_Ele12;
Int_t     trig_Ele23_Ele12;
Int_t     trig_DoubleEle33;
Int_t     trig_Mu50;
Int_t     trig_Mu17_Mu8;

Int_t aMCNLO;
Int_t nPU;
Float_t genmass;
Int_t leptType;
Float_t lept0_pt;
Float_t lept0_eta;
Float_t lept0_phi;
Float_t lept0_sceta;
Float_t lept0_miniRelIso;
Float_t lept0_trkIso;
Int_t lept0_pdgId;
Float_t lept0_normalizedChi2;
Int_t lept0_nValidMuonHits;
Int_t lept0_nMatchedStations;
Int_t lept0_nValidPixelHits;
Int_t lept0_nTrackerLayers;
Float_t lept0_muonBestTrack_dxyVTX;
Float_t lept0_muonBestTrack_dzVTX;
Float_t lept0_ptError;
Float_t lept0_sigmaIetaIeta;
Float_t lept0_dEtaIn;
Float_t lept0_dPhiIn;
Float_t lept0_hOverE;
Float_t lept0_ooEmooP;
Float_t lept0_d0;
Float_t lept0_dz;
Float_t lept0_mva;
Float_t lept0_chiso;
Float_t lept0_phoiso;
Float_t lept0_neuiso;
Float_t lept0_SIP;
Float_t lept0_sigEOverE;
Int_t lept0_expectedMissingInnerHits;
Float_t lept0_RecoSF;
Float_t lept0_SelSF;
Float_t lept0_trigSF;
Float_t lept0_SF;
Float_t lept0_RecoSF_Err;
Float_t lept0_SelSF_Err;
Float_t lept0_EnErr;
Float_t lept0_En;

Float_t lept1_pt;
Float_t lept1_eta;
Float_t lept1_phi;
Float_t lept1_sceta;
Float_t lept1_miniRelIso;
Float_t lept1_trkIso;
Int_t lept1_pdgId;
Float_t lept1_normalizedChi2;
Int_t lept1_nValidMuonHits;
Int_t lept1_nMatchedStations;
Int_t lept1_nValidPixelHits;
Int_t lept1_nTrackerLayers;
Float_t lept1_muonBestTrack_dxyVTX;
Float_t lept1_muonBestTrack_dzVTX;
Float_t lept1_ptError;
Float_t lept1_sigmaIetaIeta;
Float_t lept1_dEtaIn;
Float_t lept1_dPhiIn;
Float_t lept1_hOverE;
Float_t lept1_ooEmooP;
Float_t lept1_d0;
Float_t lept1_dz;
Float_t lept1_mva;
Float_t lept1_chiso;
Float_t lept1_phoiso;
Float_t lept1_neuiso;
Float_t lept1_SIP;
Float_t lept1_sigEOverE;
Int_t lept1_expectedMissingInnerHits;
Float_t lept1_RecoSF;
Float_t lept1_SelSF;
Float_t lept1_trigSF;
Float_t lept1_SF;
Float_t lept1_RecoSF_Err;
Float_t lept1_SelSF_Err;
Float_t lept1_EnErr;
Float_t lept1_En;
Float_t lept_dzSF;
Float_t pair_dPhi;

Float_t deltaPhi_lept;
Float_t deltaR_lept;
Int_t   isEB;
Int_t   isEE;
Float_t gamma_SF;
Float_t gamma_SF_Err;
Float_t gamma_CSEV_SF;
Float_t gamma_CSEV_SF_Err;
Float_t gamma_pt;
Float_t gamma_eta;
Float_t gamma_sceta;
Float_t gamma_phi;
Float_t gamma_HoverE;
Float_t gamma_sigmaIetaIeta;
Float_t gamma_mva;
Float_t gamma_ssmva;
Float_t gamma_R9;
Float_t gamma_ChIso;
Float_t gamma_NeuIso;
Float_t gamma_PhoIso;
Float_t gamma_ChWorstIso;
Float_t gamma_PixSeed;
Float_t gamma_EleVeto;
Float_t gamma_EnErr;
Float_t gamma_En;
Float_t deltaR_PhoLept0;
Float_t deltaR_PhoLept1;
Float_t deltaPhi_PhoLept0;
Float_t deltaPhi_PhoLept1;
Float_t gamma_geniso;
Float_t gamma_SCEtaWidth;
Float_t gamma_s4Full5x5;
Float_t gamma_sigmaIetaIphi;
Float_t gamma_SCRawE;
Float_t gamma_scetawidth;
Float_t gamma_scphiwidth;
Float_t gamma_esRR;
Float_t gamma_ESEnP1, gamma_ESEnP2;
Float_t gamma_SCEtaWidth_rw;
Float_t gamma_s4Full5x5_rw;
Float_t gamma_R9_rw;
Float_t gamma_sigmaIetaIeta_rw;
Float_t gamma_sigmaIetaIphi_rw;
Float_t etawei;

Float_t z_pt;
Float_t z_eta;
Float_t z_phi;
Float_t z_y;
Float_t z_mass;
Int_t   z_charge;
Float_t boss_pt;
Float_t boss_eta;
Float_t boss_phi;
Float_t boss_mass;
Float_t boss_massrel;



Int_t ngenEle;
Int_t ngenMu;
Int_t ngenlep;
//Int_t genlepType;
Float_t   genMuPt[maxPar];
Float_t   genMuPhi[maxPar];
Float_t   genMuEta[maxPar];
Int_t     genMuCh[maxPar];
Float_t   genlepPt[maxPar];
Float_t   genlepPhi[maxPar];
Float_t   genlepEta[maxPar];
Int_t     genlepCh[maxPar];
Float_t   genPhoEt;
Float_t   genPhoPhi;
Float_t   genPhoEta;
Float_t   genPhoY;
Float_t   gendRPhoLep1;
Float_t   gendRPhoLep2;
Float_t   gendRPhoMu1;
Float_t   gendRPhoMu2;
Float_t   genZm;
Float_t   genZpt;
Float_t   genZy;
Int_t     genZch;
Int_t     genZtype;
Float_t   genXm;
Float_t   genXpt;
Int_t     lep1fsr;
Int_t     lep2fsr;
Float_t   genCalIso03;
Float_t   genCalIso04;
Float_t   mcZm;
Float_t   lhePho1;
Float_t   lheEle1;
Float_t   lheEle2;
Float_t   lheMu1;
Float_t   lheMu2;

Float_t   lhePhoEta1;
Float_t   lheEleEta1;
Float_t   lheEleEta2;
Float_t   lheMuEta1;
Float_t   lheMuEta2;

Float_t   dRlhePhoEle1;
Float_t   dRlhePhoEle2;
Float_t   dRlhePhoMu1;
Float_t   dRlhePhoMu2;

Int_t     no_fsr_pho;
Float_t   fsr_gen_pho_pt[50];
Float_t   fsr_pho_lep_dr1[50];
Float_t   fsr_pho_lep_dr2[50];




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
  tree->Branch("rho", &rho_);

  tree->Branch("trig_doublePho60",        &trig_doublePho60);
  tree->Branch("trig_TriPho303010",       &trig_TriPho303010);
  tree->Branch("trig_TriPho35355",        &trig_TriPho35355);
  tree->Branch("trig_Ele27_WPTight",      &trig_Ele27_WPTight);
  tree->Branch("trig_Ele17_Ele12",        &trig_Ele17_Ele12);
  tree->Branch("trig_Ele23_Ele12",        &trig_Ele23_Ele12);
  tree->Branch("trig_DoubleEle33",        &trig_DoubleEle33);
  tree->Branch("trig_Mu50",               &trig_Mu50);
  tree->Branch("trig_Mu17_Mu8",           &trig_Mu17_Mu8);

  tree->Branch("leptType",                &leptType);
  tree->Branch("genmass",                 &genmass);
  tree->Branch("lept0_pt",                &lept0_pt);
  tree->Branch("lept0_eta",               &lept0_eta);
  tree->Branch("lept0_phi",               &lept0_phi);
  tree->Branch("lept0_sceta",             &lept0_sceta);
  tree->Branch("lept0_miniRelIso",        &lept0_miniRelIso);
  tree->Branch("lept0_trkIso",            &lept0_trkIso);
  tree->Branch("lept0_pdgId",             &lept0_pdgId);
  tree->Branch("lept0_normalizedChi2",    &lept0_normalizedChi2);
  tree->Branch("lept0_nValidMuonHits",    &lept0_nValidMuonHits);
  tree->Branch("lept0_nMatchedStations",  &lept0_nMatchedStations);
  tree->Branch("lept0_nValidPixelHits",   &lept0_nValidPixelHits);
  tree->Branch("lept0_nTrackerLayers",    &lept0_nTrackerLayers);
  tree->Branch("lept0_muonBestTrack_dxyVTX", &lept0_muonBestTrack_dxyVTX);
  tree->Branch("lept0_muonBestTrack_dzVTX",  &lept0_muonBestTrack_dzVTX);
  tree->Branch("lept0_ptError",              &lept0_ptError);
  tree->Branch("lept0_sigmaIetaIeta",        &lept0_sigmaIetaIeta);
  tree->Branch("lept0_dEtaIn",               &lept0_dEtaIn);
  tree->Branch("lept0_dPhiIn",               &lept0_dPhiIn);
  tree->Branch("lept0_hOverE",               &lept0_hOverE);
  tree->Branch("lept0_ooEmooP",              &lept0_ooEmooP);
  tree->Branch("lept0_d0",                   &lept0_d0);
  tree->Branch("lept0_dz",                   &lept0_dz);
  tree->Branch("lept0_mva",                  &lept0_mva);
  tree->Branch("lept0_chiso",                &lept0_chiso);
  tree->Branch("lept0_phoiso",               &lept0_phoiso);
  tree->Branch("lept0_neuiso",               &lept0_neuiso);
  tree->Branch("lept0_SIP",                  &lept0_SIP);
  tree->Branch("lept0_sigEOverE",            &lept0_sigEOverE);
  tree->Branch("lept0_expectedMissingInnerHits", &lept0_expectedMissingInnerHits);
  tree->Branch("lept0_RecoSF",               &lept0_RecoSF);
  tree->Branch("lept0_SelSF",                &lept0_SelSF);
  tree->Branch("lept0_EnErr",                &lept0_EnErr);
  tree->Branch("lept0_En",                   &lept0_En);
  tree->Branch("lept0_trigSF",               &lept0_trigSF);
  tree->Branch("lept0_SF",                   &lept0_SF);
  tree->Branch("lept0_RecoSF_Err",           &lept0_RecoSF_Err);
  tree->Branch("lept0_SelSF_Err",            &lept0_SelSF_Err);

  tree->Branch("lept1_pt",                &lept1_pt);
  tree->Branch("lept1_eta",               &lept1_eta);
  tree->Branch("lept1_phi",               &lept1_phi);
  tree->Branch("lept1_sceta",             &lept1_sceta);
  tree->Branch("lept1_miniRelIso",        &lept1_miniRelIso);
  tree->Branch("lept1_trkIso",            &lept1_trkIso);
  tree->Branch("lept1_pdgId",             &lept1_pdgId);
  tree->Branch("lept1_normalizedChi2",    &lept1_normalizedChi2);
  tree->Branch("lept1_nValidMuonHits",    &lept1_nValidMuonHits);
  tree->Branch("lept1_nMatchedStations",  &lept1_nMatchedStations);
  tree->Branch("lept1_nValidPixelHits",   &lept1_nValidPixelHits);
  tree->Branch("lept1_nTrackerLayers",    &lept1_nTrackerLayers);
  tree->Branch("lept1_muonBestTrack_dxyVTX", &lept1_muonBestTrack_dxyVTX);
  tree->Branch("lept1_muonBestTrack_dzVTX",  &lept1_muonBestTrack_dzVTX);
  tree->Branch("lept1_ptError",              &lept1_ptError);
  tree->Branch("lept1_sigmaIetaIeta",        &lept1_sigmaIetaIeta);
  tree->Branch("lept1_dEtaIn",               &lept1_dEtaIn);
  tree->Branch("lept1_dPhiIn",               &lept1_dPhiIn);
  tree->Branch("lept1_hOverE",               &lept1_hOverE);
  tree->Branch("lept1_ooEmooP",              &lept1_ooEmooP);
  tree->Branch("lept1_d0",                   &lept1_d0);
  tree->Branch("lept1_dz",                   &lept1_dz);
  tree->Branch("lept1_mva",                  &lept1_mva);
  tree->Branch("lept1_chiso",                &lept1_chiso);
  tree->Branch("lept1_phoiso",               &lept1_phoiso);
  tree->Branch("lept1_neuiso",               &lept1_neuiso);
  tree->Branch("lept1_SIP",                  &lept1_SIP);
  tree->Branch("lept1_sigEOverE",            &lept1_sigEOverE);
  tree->Branch("lept1_expectedMissingInnerHits", &lept1_expectedMissingInnerHits);
  tree->Branch("lept1_RecoSF",               &lept1_RecoSF);
  tree->Branch("lept1_SelSF",                &lept1_SelSF);
  tree->Branch("lept1_EnErr",                &lept1_EnErr);
  tree->Branch("lept1_En",                   &lept1_En);
  tree->Branch("lept1_trigSF",               &lept1_trigSF);
  tree->Branch("lept1_SF",                   &lept1_SF);
  tree->Branch("lept1_RecoSF_Err",           &lept1_RecoSF_Err);
  tree->Branch("lept1_SelSF_Err",            &lept1_SelSF_Err);
  tree->Branch("lept_dzSF",                  &lept_dzSF);
  tree->Branch("pair_dPhi",                  &pair_dPhi);

  tree->Branch("deltaR_lept",    &deltaR_lept);
  tree->Branch("deltaPhi_lept",  &deltaPhi_lept);
  tree->Branch("isEB",           &isEB);
  tree->Branch("isEE",           &isEE);
  tree->Branch("gamma_SF",       &gamma_SF);
  tree->Branch("gamma_SF_Err",   &gamma_SF_Err);
  tree->Branch("gamma_CSEV_SF",       &gamma_CSEV_SF);
  tree->Branch("gamma_CSEV_SF_Err",   &gamma_CSEV_SF_Err);

  tree->Branch("gamma_pt",       &gamma_pt);
  tree->Branch("gamma_eta",      &gamma_eta);
  tree->Branch("gamma_sceta",    &gamma_sceta);
  tree->Branch("gamma_phi",      &gamma_phi);
  tree->Branch("gamma_HoverE",   &gamma_HoverE);
  tree->Branch("gamma_sigmaIetaIeta", &gamma_sigmaIetaIeta);
  tree->Branch("gamma_mva",           &gamma_mva);
  tree->Branch("gamma_ssmva",         &gamma_ssmva);
  tree->Branch("gamma_R9",            &gamma_R9);
  tree->Branch("gamma_ChIso",         &gamma_ChIso);
  tree->Branch("gamma_NeuIso",         &gamma_NeuIso);
  tree->Branch("gamma_PhoIso",         &gamma_PhoIso);
  tree->Branch("gamma_ChWorstIso",    &gamma_ChWorstIso);
  tree->Branch("gamma_PixSeed",       &gamma_PixSeed);
  tree->Branch("gamma_EleVeto",       &gamma_EleVeto);
  tree->Branch("gamma_EnErr",         &gamma_EnErr);
  tree->Branch("gamma_En",            &gamma_En);
  tree->Branch("gamma_geniso",        &gamma_geniso);
  tree->Branch("gamma_SCEtaWidth",    &gamma_SCEtaWidth);
  tree->Branch("gamma_s4Full5x5",     &gamma_s4Full5x5);
  tree->Branch("gamma_sigmaIetaIphi", &gamma_sigmaIetaIphi);
  tree->Branch("gamma_SCRawE",        &gamma_SCRawE);
  tree->Branch("gamma_scphiwidth",    &gamma_scphiwidth);
  tree->Branch("gamma_esRR",          &gamma_esRR);
  tree->Branch("gamma_ESEnP1",        &gamma_ESEnP1);
  tree->Branch("gamma_ESEnP2",        &gamma_ESEnP2);
  tree->Branch("gamma_SCEtaWidth_rw", &gamma_SCEtaWidth_rw);
  tree->Branch("gamma_s4Full5x5_rw",  &gamma_s4Full5x5_rw);
  tree->Branch("gamma_R9_rw",         &gamma_R9_rw);
  tree->Branch("gamma_sigmaIetaIeta_rw", &gamma_sigmaIetaIeta_rw);
  tree->Branch("gamma_sigmaIetaIphi_rw", &gamma_sigmaIetaIphi_rw);
  tree->Branch("etawei",                &etawei);
  
  tree->Branch("deltaR_PhoLept0",     &deltaR_PhoLept0);
  tree->Branch("deltaR_PhoLept1",     &deltaR_PhoLept1);
  tree->Branch("deltaPhi_PhoLept0",   &deltaPhi_PhoLept0);
  tree->Branch("deltaPhi_PhoLept1",   &deltaPhi_PhoLept1);
  tree->Branch("z_pt",                &z_pt);
  tree->Branch("z_eta",               &z_eta);
  tree->Branch("z_phi",               &z_phi);
  tree->Branch("z_y",                 &z_y);
  tree->Branch("z_mass",              &z_mass);
  tree->Branch("z_charge",            &z_charge);
  tree->Branch("boss_pt",             &boss_pt);
  tree->Branch("boss_eta",            &boss_eta);
  tree->Branch("boss_phi",            &boss_phi);
  tree->Branch("boss_mass",           &boss_mass);
  tree->Branch("boss_massrel",        &boss_massrel);


  tree->Branch("ngenEle",     &ngenEle);
  tree->Branch("ngenMu",      &ngenMu);
  tree->Branch("ngenlep",      &ngenlep);
  //tree->Branch("genlepType",   &genlepType);
  //tree->Branch("genMuPt",     &genMuPt,    "genMuPt[ngenMu]/F");
  //tree->Branch("genMuPhi",    &genMuPhi,   "genMuPhi[ngenMu]/F");
  //tree->Branch("genMuEta",    &genMuEta,   "genMuEta[ngenMu]/F");
  //tree->Branch("genMuCh",     &genMuCh,    "genMuCh[ngenMu]/I");
  tree->Branch("genlepPt",    &genlepPt,   "genlepPt[ngenlep]/F");
  tree->Branch("genlepPhi",   &genlepPhi,  "genlepPhi[ngenlep]/F");
  tree->Branch("genlepEta",   &genlepEta,  "genlepEta[ngenlep]/F");
  tree->Branch("genlepCh",    &genlepCh,   "genlepCh[ngenlep]/I");
  tree->Branch("genPhoEt",    &genPhoEt,   "genPhoEt/F");
  tree->Branch("genPhoPhi",   &genPhoPhi,  "genPhoPhi/F");
  tree->Branch("genPhoEta",   &genPhoEta,  "genPhoEta/F");
  tree->Branch("genPhoY",     &genPhoY,    "genPhoY/F");
  tree->Branch("gendRPhoLep1",&gendRPhoLep1);
  tree->Branch("gendRPhoLep2",&gendRPhoLep2);
  tree->Branch("gendRPhoMu1", &gendRPhoMu1);
  tree->Branch("gendRPhoMu2", &gendRPhoMu2);
  tree->Branch("genZm",       &genZm);
  tree->Branch("genZpt",      &genZpt);
  tree->Branch("genZy",       &genZy);
  tree->Branch("genZch",      &genZch);
  tree->Branch("genZtype",    &genZtype);
  tree->Branch("genXm",       &genXm);
  tree->Branch("genXpt",      &genXpt);
  tree->Branch("genCalIso03", &genCalIso03);
  tree->Branch("genCalIso04", &genCalIso04);


}

void inittreeGen(TTree *tree){

  tree->Branch("run",         &run_);
  tree->Branch("lumi",        &lumi_);
  tree->Branch("event",       &event_);
  tree->Branch("genWeight",   &genWeight_);
  tree->Branch("puweigj",     &puweigj);
  tree->Branch("puweigj_69p2nb",     &puweigj_69p2nb);
  tree->Branch("puweigj_65nb",       &puweigj_65nb);
  tree->Branch("puweigj_63nb",       &puweigj_63nb);

  tree->Branch("ngenEle",     &ngenEle);
  tree->Branch("ngenMu",      &ngenMu);
  tree->Branch("ngenlep",      &ngenlep);
  //tree->Branch("genlepType",   &genlepType);
  //tree->Branch("genMuPt",     &genMuPt,    "genMuPt[ngenMu]/F");
  //tree->Branch("genMuPhi",    &genMuPhi,   "genMuPhi[ngenMu]/F");
  //tree->Branch("genMuEta",    &genMuEta,   "genMuEta[ngenMu]/F");
  //tree->Branch("genMuCh",     &genMuCh,    "genMuCh[ngenMu]/I");
  tree->Branch("genlepPt",    &genlepPt,   "genlepPt[ngenlep]/F");
  tree->Branch("genlepPhi",   &genlepPhi,  "genlepPhi[ngenlep]/F");
  tree->Branch("genlepEta",   &genlepEta,  "genlepEta[ngenlep]/F");
  tree->Branch("genlepCh",    &genlepCh,   "genlepCh[ngenlep]/I");
  tree->Branch("genPhoEt",    &genPhoEt,   "genPhoEt/F");
  tree->Branch("genPhoPhi",   &genPhoPhi,  "genPhoPhi/F");
  tree->Branch("genPhoEta",   &genPhoEta,  "genPhoEta/F");
  tree->Branch("gendRPhoLep1",&gendRPhoLep1);
  tree->Branch("gendRPhoLep2",&gendRPhoLep2);
  tree->Branch("gendRPhoMu1", &gendRPhoMu1);
  tree->Branch("gendRPhoMu2", &gendRPhoMu2);
  tree->Branch("genZm",       &genZm);
  tree->Branch("genZpt",      &genZpt);
  tree->Branch("genZy",       &genZy);
  tree->Branch("genZtype",    &genZtype);
  tree->Branch("genXm",       &genXm);  
  tree->Branch("genXpt",      &genXpt);
  tree->Branch("lep1fsr",     &lep1fsr);
  tree->Branch("lep2fsr",     &lep2fsr);
  tree->Branch("genCalIso03", &genCalIso03);
  tree->Branch("genCalIso04", &genCalIso04);
  tree->Branch("mcZm",        &mcZm);
  tree->Branch("lhePho1",     &lhePho1);
  tree->Branch("lheEle1",     &lheEle1);
  tree->Branch("lheEle2",     &lheEle2);
  tree->Branch("lheMu1",      &lheMu1);
  tree->Branch("lheMu2",      &lheMu2);

  tree->Branch("lhePhoEta1",    &lhePhoEta1);
  tree->Branch("lheEleEta1",    &lheEleEta1);
  tree->Branch("lheEleEta2",    &lheEleEta2);
  tree->Branch("lheMuEta1",     &lheMuEta1);
  tree->Branch("lheMuEta2",     &lheMuEta2);

  tree->Branch("dRlhePhoEle1",  &dRlhePhoEle1);
  tree->Branch("dRlhePhoEle2",  &dRlhePhoEle2);
  tree->Branch("dRlhePhoMu1",   &dRlhePhoMu1);
  tree->Branch("dRlhePhoMu2",   &dRlhePhoMu2);

  tree->Branch("no_fsr_pho",  &no_fsr_pho);
  tree->Branch("fsr_gen_pho_pt", &fsr_gen_pho_pt, "fsr_gen_pho_pt[no_fsr_pho]/F");
  tree->Branch("fsr_pho_lep_dr1", &fsr_pho_lep_dr1, "fsr_pho_lep_dr1[no_fsr_pho]/F");
  tree->Branch("fsr_pho_lep_dr2", &fsr_pho_lep_dr2, "fsr_pho_lep_dr2[no_fsr_pho]/F");
  


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
Float_t* mcMomPhi;
Float_t* mcMomMass;
UShort_t* mcStatusFlag;
Float_t* mcCalIsoDR03;
Float_t* mcCalIsoDR04;
Float_t genPho1;
Float_t genEle1;
Float_t genEle2;
Float_t genMu1;
Float_t genMu2;
Float_t genEleEta1;
Float_t genEleEta2;
Float_t genMuEta1;
Float_t genMuEta2;
Float_t genPhoEta1;
Float_t genElePhi1;
Float_t genElePhi2;
Float_t genMuPhi1;
Float_t genMuPhi2;
Float_t genPhoPhi1;




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
  mcMomPhi = data.GetPtrFloat("mcMomPhi");
  mcMomMass = data.GetPtrFloat("mcMomMass");
  mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  genWeight = data.GetFloat("genWeight");
  mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");
  mcCalIsoDR04 = data.GetPtrFloat("mcCalIsoDR04");

  genPho1 = data.GetFloat("genPho1");
  genEle1 = data.GetFloat("genEle1");
  genEle2 = data.GetFloat("genEle2");
  genMu1  = data.GetFloat("genMu1");
  genMu2  = data.GetFloat("genMu2");

  genPhoEta1 = data.GetFloat("genPhoEta1");
  genEleEta1 = data.GetFloat("genEleEta1");
  genEleEta2 = data.GetFloat("genEleEta2");
  genMuEta1  = data.GetFloat("genMuEta1");
  genMuEta2  = data.GetFloat("genMuEta2");

  genPhoPhi1 = data.GetFloat("genPhoPhi1");
  genElePhi1 = data.GetFloat("genElePhi1");
  genElePhi2 = data.GetFloat("genElePhi2");
  genMuPhi1  = data.GetFloat("genMuPhi1");
  genMuPhi2  = data.GetFloat("genMuPhi2");


  }

}


#endif
