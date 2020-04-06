#include "TH1.h"
#include "TH2.h"

TFile *feleRecoSF = new TFile("external/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root", "READ");
TFile *feleSF = new TFile("external/2016LegacyReReco_ElectronTight_Fall17V2.root", "READ");
TFile *fphoSF = new TFile("external/phoSelSF_Zg.root", "READ");
TFile *fphoCSEVSF = new TFile("external/CSEV_Eff_SF_2D.root", "READ");

TFile *feletrig_leg1 = new TFile("external/EleTrigger_Leg1_Ele23.root", "READ");
TFile *feletrig_leg2 = new TFile("external/EleTrigger_Leg2_Ele12.root", "READ");
TFile *feletrig_dz = new TFile("external/DZEff_SF_HLT_Ele23_Ele12_vsAbsEta.root", "READ");

TFile *fmutrig_BF_leg1 = new TFile("external/Eff_SF_vsPt_Mu17leg_RunBtoF.root", "READ");
TFile *fmutrig_BF_leg2 = new TFile("external/Eff_SF_vsPt_Mu8leg_RunBtoF.root", "READ");

TFile *fmutrig_GH_leg1 = new TFile("external/Eff_SF_vsPt_Mu17leg_RunGtoH.root", "READ");
TFile *fmutrig_GH_leg2 = new TFile("external/Eff_SF_vsPt_Mu8leg_RunGtoH.root", "READ");

TFile *fmutrig_dz = new TFile("external/DZEff_SF_HLT_Mu17_Mu8_vsAbsEta.root", "READ");

TFile *fmuSF_BF = new TFile("external/RunBCDEF_SF_ID_Muon.root", "READ");
TFile *fmuSF_GH = new TFile("external/RunGH_SF_ID_Muon.root", "READ");

TFile *fmuSFIso_BF = new TFile("external/RunBCDEF_SF_ISO_Muon.root", "READ");
TFile *fmuSFIso_GH = new TFile("external/RunGH_SF_ISO_Muon.root", "READ");

//get histogram of SFs from files 
TH2F *heleRecoSF = (TH2F*) feleRecoSF->Get("EGamma_SF2D");
TH2F *heleSelSF = (TH2F*) feleSF->Get("EGamma_SF2D");
TH2F *hphoSF = (TH2F*) fphoSF->Get("EGamma_SF2D");
TH2D *hphoCSEV = (TH2D*) fphoCSEVSF->Get("scale_factor");

TH2F *heletrigSF_leg1 = (TH2F*) feletrig_leg1->Get("EGamma_SF2D");
TH2F *heletrigSF_leg2 = (TH2F*) feletrig_leg2->Get("EGamma_SF2D");
TH2D *heletrigSF_dz = (TH2D*) feletrig_dz->Get("SF_DZ");

//for muon trigger
TH2D *hmutrigSF_BF_leg1 = (TH2D*) fmutrig_BF_leg1->Get("scale_factor");
TH2D *hmutrigSF_BF_leg2 = (TH2D*) fmutrig_BF_leg2->Get("scale_factor");
TH2D *hmutrigSF_GH_leg1 = (TH2D*) fmutrig_GH_leg1->Get("scale_factor");
TH2D *hmutrigSF_GH_leg2 = (TH2D*) fmutrig_GH_leg2->Get("scale_factor");

TH2D *hmutrigSF_dz = (TH2D*) fmutrig_dz->Get("SF_DZ");

TH2F *hmuSF_BF = (TH2F*) fmuSF_BF->Get("NUM_TightID_DEN_genTracks_eta_pt");
TH2F *hmuSF_GH = (TH2F*) fmuSF_GH->Get("NUM_TightID_DEN_genTracks_eta_pt");

TH2F *hmuSFiso_BF = (TH2F*) fmuSFIso_BF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");
TH2F *hmuSFiso_GH = (TH2F*) fmuSFIso_GH->Get("NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt");


float ele_RecoSF(float pt = 20., float sceta = 0.5) {

  int NBinsY = heleRecoSF->GetNbinsY();
  float maxPt = heleRecoSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 20.) pt = 20.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = heleRecoSF->GetYaxis()->FindBin(pt);
  binx = heleRecoSF->GetXaxis()->FindBin(sceta);

  float eleSF = heleRecoSF->GetBinContent(binx, biny);

  return eleSF;
}

float ele_RecoSF_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = heleRecoSF->GetNbinsY();
  float maxPt = heleRecoSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 25.) pt = 20.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = heleRecoSF->GetYaxis()->FindBin(pt);
  binx = heleRecoSF->GetXaxis()->FindBin(sceta);

  float SF_err = heleRecoSF->GetBinError(binx, biny);

  return SF_err;
}

float ele_SelSF(float pt = 20., float sceta = 0.5) {

  int NBinsY = heleSelSF->GetNbinsY();
  float maxPt = heleSelSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 20.) pt = 20.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heleSelSF->GetYaxis()->FindBin(pt);
  binx = heleSelSF->GetXaxis()->FindBin(sceta);

  float eleSF = heleSelSF->GetBinContent(binx, biny);

  return eleSF;
}

float ele_SelSF_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = heleSelSF->GetNbinsY();
  float maxPt = heleSelSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 20.) pt = 20.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heleSelSF->GetYaxis()->FindBin(pt);
  binx = heleSelSF->GetXaxis()->FindBin(sceta);

  float SF_err = heleSelSF->GetBinError(binx, biny);

  return SF_err;
}


float ele_TrigSF_leg1(float pt = 20., float sceta = 0.5) {

  int NBinsY = heletrigSF_leg1->GetNbinsY();
  float maxPt = heletrigSF_leg1->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 25.) pt = 25.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heletrigSF_leg1->GetYaxis()->FindBin(pt);
  binx = heletrigSF_leg1->GetXaxis()->FindBin(sceta);

  float eleSF = heletrigSF_leg1->GetBinContent(binx, biny);

  return eleSF;
}

float ele_TrigSF_leg1_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = heletrigSF_leg1->GetNbinsY();
  float maxPt = heletrigSF_leg1->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 25.) pt = 25.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heletrigSF_leg1->GetYaxis()->FindBin(pt);
  binx = heletrigSF_leg1->GetXaxis()->FindBin(sceta);

  float SF_err = heletrigSF_leg1->GetBinError(binx, biny);

  return SF_err;
}


float ele_TrigSF_leg2(float pt = 20., float sceta = 0.5) {

  int NBinsY = heletrigSF_leg2->GetNbinsY();
  float maxPt = heletrigSF_leg2->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 20.) pt = 20.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heletrigSF_leg2->GetYaxis()->FindBin(pt);
  binx = heletrigSF_leg2->GetXaxis()->FindBin(sceta);

  float eleSF = heletrigSF_leg2->GetBinContent(binx, biny);

  return eleSF;
}

float ele_TrigSF_leg2_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = heletrigSF_leg2->GetNbinsY();
  float maxPt = heletrigSF_leg2->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 20.) pt = 20.;
  if (pt > maxPt) pt = maxPt;
  int binx, biny;
  biny = heletrigSF_leg2->GetYaxis()->FindBin(pt);
  binx = heletrigSF_leg2->GetXaxis()->FindBin(sceta);

  float SF_err = heletrigSF_leg2->GetBinError(binx, biny);

  return SF_err;
}



float ele_trigSF_DZ(float eta1 = 0., float eta2 = 0.) {

  int binx, biny;
  biny = heletrigSF_dz->GetYaxis()->FindBin(fabs(eta1));
  binx = heletrigSF_dz->GetXaxis()->FindBin(fabs(eta2));

  float SF_DZ = heletrigSF_dz->GetBinContent(binx, biny);

  if (SF_DZ < 0.75) SF_DZ = heletrigSF_dz->GetBinContent(biny, binx);
  
  return (0.75+0.25*SF_DZ);

}


float ele_trigSF_DZ_err(float eta1 = 0., float eta2 = 0.) {

  int binx, biny;
  biny = heletrigSF_dz->GetYaxis()->FindBin(fabs(eta1));
  binx = heletrigSF_dz->GetXaxis()->FindBin(fabs(eta2));

  float SF_DZ = heletrigSF_dz->GetBinContent(binx, biny);
  float SF_err = heletrigSF_dz->GetBinError(binx, biny);

  if ( SF_DZ < 0.75) SF_err = heletrigSF_dz->GetBinError(biny, binx);

  //return sqrt(pow(0.75,2) + pow(0.25*SF_err,2));
  return SF_err;


}

/////////// muon part 
float mu_trigSF_DZ(float eta1 = 0., float eta2 = 0.) {

  int binx, biny;
  biny = hmutrigSF_dz->GetYaxis()->FindBin(fabs(eta1));
  binx = hmutrigSF_dz->GetXaxis()->FindBin(fabs(eta2));

  float SF_DZ = hmutrigSF_dz->GetBinContent(binx, biny);
  if (SF_DZ < 0.75) SF_DZ = heletrigSF_dz->GetBinContent(biny, binx);
  
  return (0.75+0.25*SF_DZ);

}


float mu_trigSF_DZ_err(float eta1 = 0., float eta2 = 0.) {

  int binx, biny;
  biny = hmutrigSF_dz->GetYaxis()->FindBin(fabs(eta1));
  binx = hmutrigSF_dz->GetXaxis()->FindBin(fabs(eta2));

  float SF_DZ = hmutrigSF_dz->GetBinContent(binx, biny);
  float SF_err = hmutrigSF_dz->GetBinError(binx, biny);

  if (SF_DZ < 0.75) SF_err = heletrigSF_dz->GetBinError(biny, binx);
  //return sqrt(pow(0.75,2) + pow(0.25*SF_err,2));
  return SF_err;

}



float muSF(float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  //for id
  int NBinsY_BF = hmuSF_BF->GetNbinsY();
  float maxPt_BF = hmuSF_BF->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmuSF_BF->GetYaxis()->FindBin(pt);
  binx_BF = hmuSF_BF->GetXaxis()->FindBin(eta);
  float muSF_BF = hmuSF_BF->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmuSF_BF->GetBinError(binx_BF, biny_BF);

  int NBinsY_GH = hmuSF_GH->GetNbinsY();
  float maxPt_GH = hmuSF_GH->GetYaxis()->GetBinCenter(NBinsY_GH);
  biny_GH = hmuSF_GH->GetYaxis()->FindBin(pt);
  binx_GH = hmuSF_GH->GetXaxis()->FindBin(eta);
  float muSF_GH = hmuSF_GH->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmuSF_GH->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  //for iso
  biny_BF = hmuSFiso_BF->GetYaxis()->FindBin(pt);
  binx_BF = hmuSFiso_BF->GetXaxis()->FindBin(eta);
  float muSFiso_BF = hmuSFiso_BF->GetBinContent(binx_BF, biny_BF);
  float Err_muSFiso_BF = hmuSFiso_BF->GetBinError(binx_BF, biny_BF);

  biny_GH = hmuSFiso_GH->GetYaxis()->FindBin(pt);
  binx_GH = hmuSFiso_GH->GetXaxis()->FindBin(eta);
  float muSFiso_GH = hmuSFiso_GH->GetBinContent(binx_GH, biny_GH);
  float Err_muSFiso_GH = hmuSFiso_GH->GetBinError(binx_GH, biny_GH);

  float muSFiso = 0.55*muSFiso_BF+ 0.45*muSFiso_GH; //lumi weight
  float Err_muSFiso = sqrt( pow(0.55*Err_muSFiso_BF,2) + pow(0.45*Err_muSFiso_GH,2)); //lumi weight

  return muSF*muSFiso;
}

float muSF_err(float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  int NBinsY_BF = hmuSF_BF->GetNbinsY();
  float maxPt_BF = hmuSF_BF->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmuSF_BF->GetYaxis()->FindBin(pt);
  binx_BF = hmuSF_BF->GetXaxis()->FindBin(eta);
  float muSF_BF = hmuSF_BF->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmuSF_BF->GetBinError(binx_BF, biny_BF);

  int NBinsY_GH = hmuSF_GH->GetNbinsY();
  float maxPt_GH = hmuSF_GH->GetYaxis()->GetBinCenter(NBinsY_GH);
  biny_GH = hmuSF_GH->GetYaxis()->FindBin(pt);
  binx_GH = hmuSF_GH->GetXaxis()->FindBin(eta);
  float muSF_GH = hmuSF_GH->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmuSF_GH->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  //for iso
  biny_BF = hmuSFiso_BF->GetYaxis()->FindBin(pt);
  binx_BF = hmuSFiso_BF->GetXaxis()->FindBin(eta);
  float muSFiso_BF = hmuSFiso_BF->GetBinContent(binx_BF, biny_BF);
  float Err_muSFiso_BF = hmuSFiso_BF->GetBinError(binx_BF, biny_BF);

  biny_GH = hmuSFiso_GH->GetYaxis()->FindBin(pt);
  binx_GH = hmuSFiso_GH->GetXaxis()->FindBin(eta);
  float muSFiso_GH = hmuSFiso_GH->GetBinContent(binx_GH, biny_GH);
  float Err_muSFiso_GH = hmuSFiso_GH->GetBinError(binx_GH, biny_GH);

  float muSFiso = 0.55*muSFiso_BF+ 0.45*muSFiso_GH; //lumi weight
  float Err_muSFiso = sqrt( pow(0.55*Err_muSFiso_BF,2) + pow(0.45*Err_muSFiso_GH,2)); //lumi weight

  float toterr = muSF*muSFiso*sqrt(pow(Err_muSF/muSF,2)+pow(Err_muSFiso/muSFiso,2));
  
  return toterr;
}


float muTrigSF_leg1 ( float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  int NBinsY_BF = hmutrigSF_BF_leg1->GetNbinsY();
  float maxPt_BF = hmutrigSF_BF_leg1->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmutrigSF_BF_leg1->GetYaxis()->FindBin(pt);
  binx_BF = hmutrigSF_BF_leg1->GetXaxis()->FindBin(fabs(eta));
  float muSF_BF = hmutrigSF_BF_leg1->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmutrigSF_BF_leg1->GetBinError(binx_BF, biny_BF);


  int NBinsY_GH = hmutrigSF_GH_leg1->GetNbinsY();
  float maxPt_GH = hmutrigSF_GH_leg1->GetYaxis()->GetBinCenter(NBinsY_GH);
  //if ( pt > maxPt_GH) pt = maxPt_GH;
  //if ( pt < 20.) pt = 20.;
  biny_GH = hmutrigSF_GH_leg1->GetYaxis()->FindBin(pt);
  binx_GH = hmutrigSF_GH_leg1->GetXaxis()->FindBin(fabs(eta));
  float muSF_GH = hmutrigSF_GH_leg1->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmutrigSF_GH_leg1->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  return muSF;

}

float muTrigSF_leg1_err ( float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  int NBinsY_BF = hmutrigSF_BF_leg1->GetNbinsY();
  float maxPt_BF = hmutrigSF_BF_leg1->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmutrigSF_BF_leg1->GetYaxis()->FindBin(pt);
  binx_BF = hmutrigSF_BF_leg1->GetXaxis()->FindBin(fabs(eta));
  float muSF_BF = hmutrigSF_BF_leg1->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmutrigSF_BF_leg1->GetBinError(binx_BF, biny_BF);


  int NBinsY_GH = hmutrigSF_GH_leg1->GetNbinsY();
  float maxPt_GH = hmutrigSF_GH_leg1->GetYaxis()->GetBinCenter(NBinsY_GH);
  //if ( pt > maxPt_GH) pt = maxPt_GH;
  //if ( pt < 20.) pt = 20.;
  biny_GH = hmutrigSF_GH_leg1->GetYaxis()->FindBin(pt);
  binx_GH = hmutrigSF_GH_leg1->GetXaxis()->FindBin(fabs(eta));
  float muSF_GH = hmutrigSF_GH_leg1->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmutrigSF_GH_leg1->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  return Err_muSF;

}


float muTrigSF_leg2 ( float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  int NBinsY_BF = hmutrigSF_BF_leg2->GetNbinsY();
  float maxPt_BF = hmutrigSF_BF_leg2->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmutrigSF_BF_leg2->GetYaxis()->FindBin(pt);
  binx_BF = hmutrigSF_BF_leg2->GetXaxis()->FindBin(fabs(eta));
  float muSF_BF = hmutrigSF_BF_leg2->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmutrigSF_BF_leg2->GetBinError(binx_BF, biny_BF);

  biny_GH = hmutrigSF_GH_leg2->GetYaxis()->FindBin(pt);
  binx_GH = hmutrigSF_GH_leg2->GetXaxis()->FindBin(fabs(eta));
  //int NBinsY_GH = hmutrigSF_GH_leg2->GetNbinsY();
  //float maxPt_GH = hmutrigSF_GH_leg2->GetYaxis()->GetBinCenter(NBinsY_GH);
  float muSF_GH = hmutrigSF_GH_leg2->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmutrigSF_GH_leg2->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  return muSF;

}

float muTrigSF_leg2_err (float pt = 20., float eta = 0.5) {

  int binx_BF, biny_BF;
  int binx_GH, biny_GH;

  int NBinsY_BF = hmutrigSF_BF_leg2->GetNbinsY();
  float maxPt_BF = hmutrigSF_BF_leg2->GetYaxis()->GetBinCenter(NBinsY_BF);
  if ( pt > maxPt_BF) pt = maxPt_BF;
  if ( pt < 20.) pt = 20.;
  biny_BF = hmutrigSF_BF_leg2->GetYaxis()->FindBin(pt);
  binx_BF = hmutrigSF_BF_leg2->GetXaxis()->FindBin(fabs(eta));
  float muSF_BF = hmutrigSF_BF_leg2->GetBinContent(binx_BF, biny_BF);
  float Err_muSF_BF = hmutrigSF_BF_leg2->GetBinError(binx_BF, biny_BF);

  biny_GH = hmutrigSF_GH_leg2->GetYaxis()->FindBin(pt);
  binx_GH = hmutrigSF_GH_leg2->GetXaxis()->FindBin(fabs(eta));
  //int NBinsY_GH = hmutrigSF_GH_leg2->GetNbinsY();
  //float maxPt_GH = hmutrigSF_GH_leg2->GetYaxis()->GetBinCenter(NBinsY_GH);
  float muSF_GH = hmutrigSF_GH_leg2->GetBinContent(binx_GH, biny_GH);
  float Err_muSF_GH = hmutrigSF_GH_leg2->GetBinError(binx_GH, biny_GH);

  float muSF = 0.55*muSF_BF+ 0.45*muSF_GH; //lumi weight 
  float Err_muSF = sqrt( pow(0.55*Err_muSF_BF,2) + pow(0.45*Err_muSF_GH,2)); //lumi weight 

  return Err_muSF;
}



//// photon part
float phoSF(float pt = 20., float sceta = 0.5) {

  int NBinsY = hphoSF->GetNbinsY();
  float maxPt = hphoSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 15.) pt = 15.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = hphoSF->GetYaxis()->FindBin(pt);
  binx = hphoSF->GetXaxis()->FindBin(sceta);

  float SF = hphoSF->GetBinContent(binx, biny);

  return SF;
}

float phoSF_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = hphoSF->GetNbinsY();
  float maxPt = hphoSF->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 15.) pt = 15.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = hphoSF->GetYaxis()->FindBin(pt);
  binx = hphoSF->GetXaxis()->FindBin(sceta);

  float SF_err = hphoSF->GetBinError(binx, biny);

  return SF_err;
}


float csev_SF(float pt = 20., float sceta = 0.5) {

  int NBinsY = hphoCSEV->GetNbinsY();
  float maxPt = hphoCSEV->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 15.) pt = 15.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = hphoCSEV->GetYaxis()->FindBin(pt);
  binx = hphoCSEV->GetXaxis()->FindBin(fabs(sceta));

  float CSEV_SF = hphoCSEV->GetBinContent(binx, biny);

  return CSEV_SF;
}

float csev_SF_err(float pt = 20., float sceta = 0.5) {

  int NBinsY = hphoCSEV->GetNbinsY();
  float maxPt = hphoCSEV->GetYaxis()->GetBinCenter(NBinsY);
  if (pt < 15.) pt = 15.;
  if (pt > maxPt) pt = maxPt;

  int binx, biny;
  biny = hphoCSEV->GetYaxis()->FindBin(pt);
  binx = hphoCSEV->GetXaxis()->FindBin(fabs(sceta));

  float CSEV_SF_err = hphoCSEV->GetBinError(binx, biny);

  return CSEV_SF_err;
}
