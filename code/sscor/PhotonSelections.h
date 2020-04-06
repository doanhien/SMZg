#include <algorithm>
#include <TMath.h>
#include <TMVA/Reader.h>

#include "TGraph.h"


bool preselect_pho_Zg(int ipho) {

  bool pass = false;
  if ( fabs(phoSCEta[ipho]) < 1.5 
       && phoHoverE[ipho] < 0.08 
       && phoSigmaIEtaIEtaFull5x5[ipho] < 0.015 //0.012
       && phoPFPhoIso[ipho] < 15.
       && phoPFChWorstIso[ipho] < 15.) 
    pass = true;
  else if ( fabs(phoSCEta[ipho]) > 1.5
	    &&phoHoverE[ipho] < 0.05
	    && phoSigmaIEtaIEtaFull5x5[ipho] < 0.045 //0.027
	    && phoPFPhoIso[ipho] < 15.
	    && phoPFChWorstIso[ipho] < 15.)
    pass = true;

  return pass;
}

/*
void preselect_pho_Zg(vector<int> &selectedPho) {

  selectedPho.clear();
  vector<int> passUnsorted;
  passUnsorted.clear();
  vector<float> phoPtUnsorted;
  phoPtUnsorted.clear();

  for (int ipho = 0; ipho < nPho; ipho++) {
  if (phoCalibEt[ipho] < 15. ) continue;
  if (fabs(phoSCEta[ipho]) > 2.5) continue;
  if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;
  //preselection
  if ( fabs(phoSCEta[ipho]) < 1.5 ) {
    if (phoHoverE[ipho] > 0.08 ) continue;
    if (phoSigmaIEtaIEtaFull5x5[ipho] > 0.015) continue;
    if (phoPFPhoIso[ipho] > 15.) continue;
    if (phoPFChWorstIso[ipho] > 15.) continue;
  }
  else if ( fabs(phoSCEta[ipho]) > 1.5 ) {
    if (phoHoverE[ipho] > 0.05 ) continue;
    if (phoSigmaIEtaIEtaFull5x5[ipho] > 0.045) continue; //0.027
    if (phoPFPhoIso[ipho] > 15.) continue;
    if (phoPFChWorstIso[ipho] > 15.) continue;
  }

  if (phoEleVeto[ipho] != 1) continue;

  passUnsorted.push_back(ipho);
  phoPtUnsorted.push_back(phoCalibEt[ipho]);
  }

  //sorted pt
  int siz = (int) passUnsorted.size();
  if (siz < 1) return;

  int ind[siz];
  TMath::Sort(siz, &phoPtUnsorted.front(), ind);
  for (int i = 0; i < siz; ++i) {
    selectedPho.push_back(passUnsorted[ind[i]]);
  }

}
*/

float PhotonSSMVA_rw(TreeReader &data, Int_t i, TGraph *tgr[20]) {

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
  //Float_t* phoESEn           = data.GetPtrFloat("phoESEn");
  Float_t* phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");
  Float_t* phoESEnP1          = data.GetPtrFloat("phoESEnP1");
  Float_t* phoESEnP2          = data.GetPtrFloat("phoESEnP2");

  Float_t* phoSigmaIEtaIEtaFull5x5  = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  Float_t* phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
  Float_t* phoE2x2Full5x5           = data.GetPtrFloat("phoE2x2Full5x5");
  Float_t* phoE5x5Full5x5           = data.GetPtrFloat("phoE5x5Full5x5");


  static float phoEt_, phoEta_, phoPhi_, phoR9_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoSCEta_, phoSCRawE_;
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
    tmvaReader[iBE]->AddVariable( "s4 := e2x2Full5x5/e5x5Full5x5",    &s4Full5x5 );
    tmvaReader[iBE]->AddVariable("recoSCEta", &phoSCEta_);
    tmvaReader[iBE]->AddVariable("rawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("scEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("scPhiWidth", &phoSCPhiWidth_);
    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("ESEn := esEn/rawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("esRR", &phoESEffSigmaRR_);
    }
    tmvaReader[iBE]->AddVariable("rho", &rho_);
    //if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/TMVAnalysis_BDT_EB_v1.weights.xml");
    //else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/TMVAnalysis_BDT_EE_v1.weights.xml");
    //if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/420/TMVAnalysis_BDT_EB.weights.xml");
    //else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/420/TMVAnalysis_BDT_EE.weights.xml");
    if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/upto6000/TMVAnalysis_BDT_EB.weights.xml");
    else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/upto6000/TMVAnalysis_BDT_EE.weights.xml");
  }

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

  /*
  //eta weighted for endcap to check MVA
  float eta_wei = 1;
  if (TMath::Abs(phoSCEta[i])>1.5) {
    if (phoSCEta[i]>-2.5 && phoSCEta[i]<-2.3)  eta_wei = 2.422;
    else if ( phoSCEta[i]>-2.3 && phoSCEta[i]<-2.1) eta_wei = 1.194;
    else if ( phoSCEta[i]>-2.1 && phoSCEta[i]<-1.9) eta_wei = 1.343;
    else if ( phoSCEta[i]>-1.9 && phoSCEta[i]<-1.7) eta_wei = 0.617;
    else if ( phoSCEta[i]>-1.7 && phoSCEta[i]<-1.5) eta_wei = 0.474;
    else if ( phoSCEta[i]> 1.5 && phoSCEta[i]< 1.7) eta_wei = 0.678;
    else if ( phoSCEta[i]> 1.7 && phoSCEta[i]< 1.9) eta_wei = 0.516;
    else if ( phoSCEta[i]> 1.9 && phoSCEta[i]< 2.1) eta_wei = 1.070;
    else if ( phoSCEta[i]> 2.1 && phoSCEta[i]< 2.3) eta_wei = 1.444;
    else if ( phoSCEta[i]> 2.1 && phoSCEta[i]< 2.5) eta_wei = 2.825;
    phoSCEta_ *= eta_wei;
  }
  */
  
  //shower shape correction for mc
    if(isData!=1){
    if(TMath::Abs(phoSCEta[i])<1.5) {
      //float etaw_rw  = tgr[3]->Eval(phoSCEtaWidth_);
      //phoSCEtaWidth_    = tgr[0]->Eval(etaw_rw);
      phoSCEtaWidth_    = tgr[0]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[1]->Eval(s4Full5x5);
      phoR9_            = tgr[2]->Eval(phoR9_);
      sieieFull5x5      = tgr[3]->Eval(sieieFull5x5);
      sieipFull5x5      = tgr[4]->Eval(sieipFull5x5);
      //rho_              = tgr[5]->Eval(rho);
      //phoSCPhiWidth_    = tgr[6]->Eval(phoSCPhiWidth_);

      
    }else{
      //float etaw_rw  = tgr[10]->Eval(phoSCEtaWidth_);
      //phoSCEtaWidth_    = tgr[7]->Eval(etaw_rw);
      phoSCEtaWidth_    = tgr[7]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[8]->Eval(s4Full5x5);
      phoR9_            = tgr[9]->Eval(phoR9_);
      sieieFull5x5      = tgr[10]->Eval(sieieFull5x5);
      sieipFull5x5      = tgr[11]->Eval(sieipFull5x5);
      ////rho_              = tgr[12]->Eval(rho);
      //phoSCPhiWidth_    = tgr[13]->Eval(phoSCPhiWidth_);

    }
    }


  return tmvaReader[iBE]->EvaluateMVA("BDT");

}


float PhotonSSMVA(TreeReader &data, Int_t i) {

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
  //Float_t* phoESEn           = data.GetPtrFloat("phoESEn");
  Float_t* phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");
  Float_t* phoESEnP1          = data.GetPtrFloat("phoESEnP1");
  Float_t* phoESEnP2          = data.GetPtrFloat("phoESEnP2");

  Float_t* phoSigmaIEtaIEtaFull5x5  = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
  Float_t* phoSigmaIEtaIPhiFull5x5  = data.GetPtrFloat("phoSigmaIEtaIPhiFull5x5");
  Float_t* phoE2x2Full5x5           = data.GetPtrFloat("phoE2x2Full5x5");
  Float_t* phoE5x5Full5x5           = data.GetPtrFloat("phoE5x5Full5x5");


  static float phoEt_, phoEta_, phoPhi_, phoR9_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoSCEta_, phoSCRawE_;
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
    tmvaReader[iBE]->AddVariable( "s4 := e2x2Full5x5/e5x5Full5x5",    &s4Full5x5 );
    tmvaReader[iBE]->AddVariable("recoSCEta", &phoSCEta_);
    tmvaReader[iBE]->AddVariable("rawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("scEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("scPhiWidth", &phoSCPhiWidth_);
    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("ESEn := esEn/rawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("esRR", &phoESEffSigmaRR_);
    }
    tmvaReader[iBE]->AddVariable("rho", &rho_);
    if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/upto6000/TMVAnalysis_BDT_EB.weights.xml");
    else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/summer16/upto6000/TMVAnalysis_BDT_EE.weights.xml");
  }

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
  float *mcCalIsoDR03 = data.GetPtrFloat("mcCalIsoDR03");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  bool match = false;
  int genIndex = -1;

  for(int imc=0; imc<nMC; imc++){
    if (mcPID[imc]!=22) continue;
    if (((mcStatusFlag[imc]>>1)&1) != 1 && ((mcStatusFlag[imc]>>0)&1) != 1) continue;
    //if ( abs(mcMomPID[imc])>16 && mcMomPID[imc] != 21 && mcMomPID[imc] != 23 && abs(mcMomPID[imc]) != 999) continue; //no need because (mcStatusFlag[imc]>>1)&1) already count
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.1 ) {
      //if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.2 && fabs(phoEt[iPho]-mcPt[imc])/mcPt[imc] < 0.2 ) {
      match = true;
      genIndex = imc;
      //break;
    }

    if (match && mcCalIsoDR03[genIndex] < 5.) {
      //if (match) {
      genmatch = 1;
      break;
    }
    else genmatch = 0;
  }

  return genmatch;
}


bool phoMatcher (TreeReader &data, int iPho) {

  bool genmatch = false;
  if ( !data.HasMC()) genmatch = false ;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  for(int imc=0; imc<nMC; imc++){
    //if (mcPID[imc]!=22 && abs(mcPID[imc])!=11) continue;
    if (abs(mcPID[imc])!=11) continue;
    if (((mcStatusFlag[imc]>>0)&1) != 1) continue;
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.3 ) {
      genmatch = true;
      break;
    }
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

  for (Int_t i=0; i<nMC; ++i) {
    if (abs(mcPID[i]) != 11) continue;
    //if ( ((mcStatusFlag[i] >> 0) & 1) == 0 && ((mcStatusFlag[i] >> 1) & 1) == 0 ) continue;
    //mcpho.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0);
    if ( deltaR(mcEta[i],mcPhi[i],phoEta[ipho],phoPhi[ipho]) < 0.2 ) {
      match = 1;
      break;
    }
  }

  return match;

}


