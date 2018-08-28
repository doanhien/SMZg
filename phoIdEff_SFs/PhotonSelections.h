#include <algorithm>
#include <TMath.h>
#include <TMVA/Reader.h>

bool select_photonXZg( int iPho) {

  bool pass = false;

  if (phoCalibEt[iPho] > 20. && fabs(phoSCEta[iPho]) < 2.5) {
    if (fabs(phoSCEta[iPho]) < 1.479
	&& phohasPixelSeed[iPho]==0
	//&& phoEleVeto[iPho] ==1
	&& phoSigmaIEtaIEtaFull5x5[iPho] < 0.0102
	&& phoHoverE[iPho] < 0.05
	&& phoPFChIso[iPho] < 2.5)
      pass = true;
    else if (fabs(phoSCEta[iPho]) > 1.57
	     && phohasPixelSeed[iPho]==0
	     //&& phoEleVeto[iPho] ==1
	     && phoSigmaIEtaIEtaFull5x5[iPho] < 0.0274
	     && phoHoverE[iPho] < 0.05 
	     && phoPFChIso[iPho] < 2.5 )
      pass = true;
  }

  return pass;

}

bool select_photon_80X (int iPho) {

  bool pass = false;

  int index;
  if (fabs(phoSCEta[iPho]) < 1.0) index = 0;
  else if (fabs(phoSCEta[iPho])>1.0 && fabs(phoSCEta[iPho])<1.479) index = 1;
  else if (fabs(phoSCEta[iPho])>1.479 && fabs(phoSCEta[iPho])<2.0) index = 2;
  else if (fabs(phoSCEta[iPho])>2.0 && fabs(phoSCEta[iPho])<2.2) index = 3;
  else if (fabs(phoSCEta[iPho])>2.2 && fabs(phoSCEta[iPho])<2.3) index = 4;
  else if (fabs(phoSCEta[iPho])>2.3 && fabs(phoSCEta[iPho])<2.4) index = 5;
  else index = 6;

  Double_t EACh[7] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
  Double_t EANeu[7] = {0.0597, 0.0807, 0.0629, 0.0197, 0.0184, 0.0284, 0.0591};
  Double_t EAPho[7] = {0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998};

  Double_t rhoChIso = TMath::Max(phoPFChIso[iPho] - rho*EACh[index], 0.0);
  Double_t rhoNeuIso = TMath::Max((phoPFNeuIso[iPho] - rho*EANeu[index]), 0.0);
  Double_t rhoPhoIso = TMath::Max(phoPFPhoIso[iPho] - rho*EAPho[index], 0.);

  if (phoCalibEt[iPho] > 20. && fabs(phoSCEta[iPho]) < 2.5) {
    if (fabs(phoSCEta[iPho]) < 1.479
	//&& phohasPixelSeed[iPho]==0
	&& phoEleVeto[iPho] == 1  
	&& phoHoverE[iPho] < 0.0597
	&& phoSigmaIEtaIEtaFull5x5[iPho] < 0.01031
	&& rhoChIso < 1.295
	&& rhoNeuIso < 10.910+0.0148*phoCalibEt[iPho]+0.000017*pow(phoCalibEt[iPho],2) 
	&& rhoPhoIso < 3.630+0.0047*phoCalibEt[iPho])
      pass = true;
    else if ( fabs(phoSCEta[iPho]) > 1.57
	      //&& phohasPixelSeed[iPho]==0
	      && phoEleVeto[iPho] == 1
	      && phoHoverE[iPho] < 0.0481
	      && phoSigmaIEtaIEtaFull5x5[iPho] < 0.03013
	      && rhoChIso < 1.01 
	      && rhoNeuIso < 5.931+0.0163*phoCalibEt[iPho]+0.000014*pow(phoCalibEt[iPho],2)
	      && rhoPhoIso < 6.641+0.0034*phoCalibEt[iPho])
      pass = true;
  }

  return pass;

}

bool preselect_pho_Zg(int ipho) {

  bool pass = false;
  if ( fabs(phoSCEta[ipho]) < 1.5 
       && phoHoverE[ipho] < 0.08 
       //&& phoSigmaIEtaIEtaFull5x5[ipho] < 0.012
       && phoSigmaIEtaIEtaFull5x5[ipho] < 0.015
       && phoPFPhoIso[ipho] < 15.
       && phoPFChWorstIso[ipho] < 15.) 
    pass = true;
  else if ( fabs(phoSCEta[ipho]) > 1.5
	    &&phoHoverE[ipho] < 0.05
	    //&& phoSigmaIEtaIEtaFull5x5[ipho] < 0.027
	    && phoSigmaIEtaIEtaFull5x5[ipho] < 0.045
	    && phoPFPhoIso[ipho] < 15.
	    && phoPFChWorstIso[ipho] < 15.)
    pass = true;

  return pass;
}

/*
bool select_pho_mvanoIso(int ipho, TGraph *tgr[20]) {

  // classification variables
  static float phoEt_, phoEta_, phoPhi_, phoR9_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoSCEta_, phoSCRawE_;
  //static float phoPFPhoIso_, phoPFChIso_, phoPFChIsoWorst_; 
  static float phoESEnToRawE_, phoESEffSigmaRR_;

  static float sieieFull5x5, sieipFull5x5, s4Full5x5; //s13Full5x5, s25Full5x5; 

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
    tmvaReader[iBE]->AddVariable( "sieieFull5x5",                     &sieieFull5x5 );
    tmvaReader[iBE]->AddVariable( "sieipFull5x5",                     &sieipFull5x5 );
    //tmvaReader[iBE]->AddVariable( "s13 := e1x3Full5x5/e5x5Full5x5",   &s13Full5x5 ); 
    tmvaReader[iBE]->AddVariable( "s4 := e2x2Full5x5/e5x5Full5x5",    &s4Full5x5 );
    //tmvaReader[iBE]->AddVariable( "s25 := e2x5Full5x5/e5x5Full5x5",   &s25Full5x5 ); 
    tmvaReader[iBE]->AddVariable("recoSCEta", &phoSCEta_);
    tmvaReader[iBE]->AddVariable("rawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("scEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("scPhiWidth", &phoSCPhiWidth_);
    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("ESEn := esEn/rawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("esRR", &phoESEffSigmaRR_);
    }
    tmvaReader[iBE]->AddVariable("rho", &rho_);

    //read weight files
    if (iBE == 0){
      tmvaReader[0]->BookMVA("BDT", "TMVA/spring16/TMVAnalysis_BDT_EB.weights.xml");
    }else{
      tmvaReader[1]->BookMVA("BDT", "TMVA/spring16/TMVAnalysis_BDT_EE.weights.xml");
    }
  }// one-time initialization 

  // set MVA variables
  phoPhi_ = phoPhi[ipho];
  phoR9_ = phoR9[ipho];
  phoSCEta_ = phoSCEta[ipho];
  phoSCRawE_ = phoSCRawE[ipho];
  phoSCEtaWidth_ = phoSCEtaWidth[ipho];
  phoSCPhiWidth_ = phoSCPhiWidth[ipho];
  rho_ = rho;
  phoESEnToRawE_ = (phoESEnP1[ipho]+phoESEnP2[ipho])/phoSCRawE[ipho];
  phoESEffSigmaRR_= phoESEffSigmaRR[ipho];
  phoEt_ = phoEt[ipho];
  phoEta_ = phoEta[ipho];
  sieieFull5x5 = phoSigmaIEtaIEtaFull5x5[ipho];
  sieipFull5x5 = phoSigmaIEtaIPhiFull5x5[ipho];
  s4Full5x5 = phoE2x2Full5x5[ipho]/phoE5x5Full5x5[ipho];

  if(isData!=1){
    if(TMath::Abs(phoSCEta[ipho])<1.5) {
      phoSCEtaWidth_    = tgr[0]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[1]->Eval(s4Full5x5);
      phoR9_            = tgr[2]->Eval(phoR9_);
    }else{
      phoSCEtaWidth_    = tgr[3]->Eval(phoSCEtaWidth_);
      s4Full5x5         = tgr[4]->Eval(s4Full5x5);
      phoR9_            = tgr[5]->Eval(phoR9_);

    }
  }
  return tmvaReader[iBE]->EvaluateMVA("BDT");

}
*/

/*
float PhotonSSMVA(Int_t i) {

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
    if (iBE == 0) tmvaReader[0]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/spring16/TMVAnalysis_BDT_EB.weights.xml");
    else tmvaReader[1]->BookMVA("BDT", "/data4/tdoan/analyzer/Zg_SM/2016data/ana/TMVA/spring16/TMVAnalysis_BDT_EE.weights.xml");
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
*/

bool phoMatcher_prompt (TreeReader &data, int iPho) {

  bool match = false;
  if ( !data.HasMC()) match= false ;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  for(int imc=0; imc<nMC; imc++){
    //if (mcPID[imc]!=22 || ( (mcStatusFlag[imc]>>0&1) != 1 && (mcStatusFlag[imc]>>1&1) != 1) ) continue;  //for Z+jets
    //if (abs(mcPID[imc])!=11 || ( (mcStatusFlag[imc]>>0&1) != 1 && (mcStatusFlag[imc]>>1&1) != 1 )) continue; //for WZ and ZZ
    if (mcPID[imc]!=22) continue;
    if ( (mcStatusFlag[imc]>>0&1) != 1 &&  (mcStatusFlag[imc]>>0&1) != 1 ) continue; //check photon internal conversion
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.1 ) match = true ;
  }

  return match;
}
 
//match to gen ele
bool phoMatcher (TreeReader &data, int iPho) {

  bool match = false;
  if ( !data.HasMC()) match= false ;

  Int_t nMC        = data.GetInt("nMC");
  Int_t *mcPID     = data.GetPtrInt("mcPID");
  Int_t *mcMomPID  = data.GetPtrInt("mcMomPID");
  float *mcMomMass = data.GetPtrFloat("mcMomMass");
  float *mcEta     = data.GetPtrFloat("mcEta");
  float *mcPhi     = data.GetPtrFloat("mcPhi");
  float *mcMass    = data.GetPtrFloat("mcMass");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");

  for(int imc=0; imc<nMC; imc++){
    //if (mcPID[imc]!=22 ) continue;
    //if ( (mcStatusFlag[imc]>>0&1) != 1 && (mcStatusFlag[imc]>>1&1) != 1 ) continue;  //for Z+jets
    //if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.1 ) match = true ;
    if ( abs(mcPID[imc])!=11 ) continue;
    if ( (mcStatusFlag[imc]>>0&1) != 1 ) continue;
    if ( deltaR(mcEta[imc],mcPhi[imc],phoEta[iPho],phoPhi[iPho]) < 0.2 ) match = true ;
  }

  return match;
}
 
