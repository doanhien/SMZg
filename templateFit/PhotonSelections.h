#include "TMVA/Reader.h"
#include "TGraph.h"

using namespace std;

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

Int_t PassPhotonPreselections(TreeReader &data, Int_t i) {

  Int_t pass = 1;

  Int_t  nPho                    = data.GetInt("nPho");
  //float* phoEt                   = data.GetPtrFloat("phoEt");
  float* phoEt                   = data.GetPtrFloat("phoCalibEt");
  float* phoSCEta                = data.GetPtrFloat("phoSCEta");
  float* phoHoverE               = data.GetPtrFloat("phoHoverE");
  float* phoSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  float* phoPFPhoIso             = data.GetPtrFloat("phoPFPhoIso");
  float* phoPFChWorstIso         = data.GetPtrFloat("phoPFChWorstIso");
  Int_t* phoEleVeto              = data.GetPtrInt("phoEleVeto");

  Float_t cut_hoe[2]       = {0.08,  0.05};
  Float_t cut_sieie[2]     = {0.015, 0.045};
  //Float_t cut_sieie[2]     = {0.012, 0.027};
  Float_t cut_phoIso[2]    = {15,    15};
  Float_t cut_worsCHIso[2] = {15,    15};

  if (phoEt[i] < 10) pass = 0;;
  if (fabs(phoSCEta[i]) > 2.5) pass = 0;;
  if (fabs(phoSCEta[i]) > 1.4442 && fabs(phoSCEta[i]) < 1.566) pass = 0;;
  
  if (phoEleVeto[i] == 0) pass = 0;;
  
  Int_t iEB = (fabs(phoSCEta[i]) < 1.4442) ? 0 : 1;

  if (phoHoverE[i]               > cut_hoe[iEB])       pass = 0;;
  if (phoSigmaIEtaIEtaFull5x5[i] > cut_sieie[iEB])     pass = 0;; 
  if (phoPFPhoIso[i]             > cut_phoIso[iEB])    pass = 0;;
  if (phoPFChWorstIso[i]         > cut_worsCHIso[iEB]) pass = 0;;
  
  return pass;
}

float PhotonSSMVA(TreeReader &data, Int_t i, TGraph *tgr[20]) {

  /* Photon identification with the Zgamma MVA. Returns the MVA evaluated value.
   *
   * Documentation:
   * https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=298231
   *
   * data = handle providing access to an input event;
   * i = index of a photon candidate to consider.
   */

  Bool_t  isData             = data.GetBool("isData");
  // load necessary tree branches
  Float_t* phoEt             = data.GetPtrFloat("phoEt");
  Float_t* phoEta            = data.GetPtrFloat("phoEta");
  Float_t* phoPhi            = data.GetPtrFloat("phoPhi");
  Float_t* phoR9             = data.GetPtrFloat("phoR9");
  Float_t* phoSCEta          = data.GetPtrFloat("phoSCEta");
  Float_t* phoSCRawE         = data.GetPtrFloat("phoSCRawE");
  Float_t* phoSCEtaWidth     = data.GetPtrFloat("phoSCEtaWidth");
  Float_t* phoSCPhiWidth     = data.GetPtrFloat("phoSCPhiWidth");
  Float_t  rho               = data.GetFloat("rho");
  /* Float_t* phoPFPhoIso       = data.GetPtrFloat("phoPFPhoIso"); */
  /* Float_t* phoPFChIso        = data.GetPtrFloat("phoPFChIso"); */
  /* Float_t* phoPFChWorstIso   = data.GetPtrFloat("phoPFChWorstIso"); */
  Float_t* phoESEnP1           = data.GetPtrFloat("phoESEnP1");
  Float_t* phoESEnP2           = data.GetPtrFloat("phoESEnP2");
  Float_t* phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");

  Float_t* phoSigmaIEtaIEtaFull5x5  = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  Float_t* phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
  /* Float_t* phoE1x3Full5x5           = data.GetPtrFloat("phoE1x3Full5x5"); */
  Float_t* phoE2x2Full5x5           = data.GetPtrFloat("phoE2x2Full5x5");
  Float_t* phoE5x5Full5x5           = data.GetPtrFloat("phoE5x5Full5x5");
  /* Float_t* phoE2x5MaxFull5x5        = data.GetPtrFloat("phoE2x5MaxFull5x5"); */

  // classification variables
  static float phoEt_, phoEta_, phoPhi_, phoR9_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoSCEta_, phoSCRawE_;
  //static float phoPFPhoIso_, phoPFChIso_, phoPFChIsoWorst_;
  static float phoESEnToRawE_, phoESEffSigmaRR_;

  static float sieieFull5x5, sieipFull5x5, s13Full5x5, s4Full5x5, s25Full5x5;

  // MVA classifiers for 0=ECAL barrel and 1=ECAL endcaps
  static TMVA::Reader* tmvaReader[2] = {NULL, NULL};

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (fabs(phoSCEta[i]) < 1.479) ? 0 : 1;

  // one-time MVA initialization
  if (!tmvaReader[iBE]) {
    tmvaReader[iBE] = new TMVA::Reader("!Color:Silent");

    // add classification variables
    tmvaReader[iBE]->AddVariable("recoPhi", &phoPhi_);
    tmvaReader[iBE]->AddVariable("r9", &phoR9_);
    tmvaReader[iBE]->AddVariable( "sieieFull5x5",             &sieieFull5x5 );     
    tmvaReader[iBE]->AddVariable( "sieipFull5x5",             &sieipFull5x5 );     
    /* tmvaReader[iBE]->AddVariable( "s13 := e1x3Full5x5/e5x5Full5x5",   &s13Full5x5 );         */
    tmvaReader[iBE]->AddVariable( "s4 := e2x2Full5x5/e5x5Full5x5",    &s4Full5x5 );       
    /* tmvaReader[iBE]->AddVariable( "s25 := e2x5Full5x5/e5x5Full5x5",   &s25Full5x5 );         */
    tmvaReader[iBE]->AddVariable("recoSCEta", &phoSCEta_);
    tmvaReader[iBE]->AddVariable("rawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("scEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("scPhiWidth", &phoSCPhiWidth_);
    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("ESEn := esEn/rawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("esRR", &phoESEffSigmaRR_);
    }
    tmvaReader[iBE]->AddVariable("rho", &rho_);

    /*tmvaReaderVariable("phoIsoRaw", &phoPFPhoIso_);*/
      /* tmvaReader[iBE]->AddVariable("chIsoRaw", &phoPFChIso_); */
      /* tmvaReader[iBE]->AddVariable("chWorstRaw", &phoPFChIsoWorst_); */

      //tmvaReader[iBE]->AddVariable("recoPt", &phoEt_);
      // FIXME: why do we need this?
      /* tmvaReader[iBE]->AddSpectator("recoPt", &phoEt_); */
      /* tmvaReader[iBE]->AddSpectator("recoEta", &phoEta_); */

      // read weight files
      //if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/spring16/TMVAnalysis_BDT_EB.weights.xml");
      //else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/spring16/TMVAnalysis_BDT_EE.weights.xml");

    //summer16
    //if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/TMVAnalysis_BDT_EB_v1.weights.xml");
    //else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/TMVAnalysis_BDT_EE_v1.weights.xml");

      if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/420/TMVAnalysis_BDT_EB.weights.xml");
      else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/420/TMVAnalysis_BDT_EE.weights.xml");
    
  } // one-time initialization

  // set MVA variables
  phoPhi_          = phoPhi[i];
  phoR9_           = phoR9[i];
  phoSCEta_        = phoSCEta[i];
  phoSCRawE_       = phoSCRawE[i];
  phoSCEtaWidth_   = phoSCEtaWidth[i];
  phoSCPhiWidth_   = phoSCPhiWidth[i];
  rho_             = rho;
  phoESEnToRawE_   = (phoESEnP1[i] + phoESEnP2[i])/phoSCRawE[i];
  phoESEffSigmaRR_ = phoESEffSigmaRR[i];
  phoEt_           = phoEt[i];
  phoEta_          = phoEta[i];
  sieieFull5x5     = phoSigmaIEtaIEtaFull5x5[i];
  sieipFull5x5     = phoSigmaIEtaIPhiFull5x5[i];
  s4Full5x5        = phoE2x2Full5x5[i]/phoE5x5Full5x5[i];
  /* s13Full5x5    = phoE1x3Full5x5[i]/phoE5x5Full5x5[i]; */
  /* phoPFPhoIso_     = phoPFPhoIso[i]; */
  /* phoPFChIso_      = phoPFChIso[i]; */
  /* phoPFChIsoWorst_ = phoPFChWorstIso[i]; */
  /* s25Full5x5    = phoE2x5MaxFull5x5[i]/phoE5x5Full5x5[i]; */

  //showershape correction for MC
  if(isData!=1){
    if(TMath::Abs(phoSCEta[i])<1.5) {
      phoSCEtaWidth_    = tgr[0]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[1]->Eval(s4Full5x5);
      phoR9_            = tgr[2]->Eval(phoR9_);
      sieieFull5x5      = tgr[3]->Eval(sieieFull5x5);
    }else{
      phoSCEtaWidth_    = tgr[4]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[5]->Eval(s4Full5x5);
      phoR9_            = tgr[6]->Eval(phoR9_);
      sieieFull5x5      = tgr[7]->Eval(sieieFull5x5);
    }
  }
  
  return tmvaReader[iBE]->EvaluateMVA("BDT");

}


bool phoMatcher_prompt (TreeReader &data, int iPho) {

  bool genmatch = false;
  if ( !data.HasMC()) genmatch = false ;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  float *mcPt      = data.GetPtrFloat("mcPt");
  float* mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  float *phoEt = data.GetPtrFloat("phoCalibEt");
  float *phoEta = data.GetPtrFloat("phoEta");
  float *phoPhi = data.GetPtrFloat("phoPhi");

  bool match = false;
  int genIndex = -1;

  for(int imc=0; imc<nMC; imc++){
    if (mcPID[imc]!=22) continue;
    if ( ((mcStatusFlag[imc]>>0)&1) != 1 &&  ((mcStatusFlag[imc]>>1)&1) != 1 ) continue; //check photon internal conversion                                             
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.2 && fabs(phoEt[iPho]-mcPt[imc])/mcPt[imc] < 0.2 ) {
      match = true;
      genIndex = imc;
    }

    if (match && mcCalIsoDR03[genIndex] < 5.) {
      genmatch = 1;
      break;
    }
    else genmatch = 0;
  }

  return genmatch;
}

bool pho_matchele (TreeReader &data, int ipho) {

  bool match = 0;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  Float_t* mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");

  float *phoEt = data.GetPtrFloat("phoCalibEt");
  float *phoEta = data.GetPtrFloat("phoEta");
  float *phoPhi = data.GetPtrFloat("phoPhi");

  for (Int_t i=0; i<nMC; ++i) {
    if (abs(mcPID[i]) != 11) continue;
    //if ( ((mcStatusFlag[i] >> 0) & 1) == 0 && ((mcStatusFlag[i] >> 1) & 1) == 0 ) continue; 
    //mcpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
    if ( deltaR(mcEta[i],mcPhi[i],phoEta[ipho],phoPhi[ipho]) < 0.7 ) {
      match = 1;
      break;
    }
  }

  return match;

}


//matching fsr
bool pho_matchfsr (TreeReader &data, int ipho) {

  bool match = 0;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  Float_t* mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");

  float *phoEt = data.GetPtrFloat("phoCalibEt");
  float *phoEta = data.GetPtrFloat("phoEta");
  float *phoPhi = data.GetPtrFloat("phoPhi");

  for (Int_t i=0; i<nMC; ++i) {
    if (abs(mcPID[i]) != 11) continue;
    if ( ((mcStatusFlag[i] >> 0) & 1) == 0 && ((mcStatusFlag[i] >> 1) & 1) == 0 ) continue;
    //mcpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
    if ( deltaR(mcEta[i],mcPhi[i],phoEta[ipho],phoPhi[ipho]) < 0.7 ) {
      match = 1;
      break;
    }
  }

  return match;

}
