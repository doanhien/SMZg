#ifndef ElectronSelection_h__
#define ElectronSelection_h__

/*
bool pass_cutbased_80X (int iEle, int iWP) {  // iWP: 0 -veto, 1-loose, 2 - medium, 3-tight

  bool pass = true;

  Float_t sIeIeCut_EB[4]     = {0.0115,  0.011,   0.00998,  0.00998};
  Float_t dEtaInCut_EB[4]    = {0.00749, 0.00477, 0.00311,  0.00308 };
  Float_t dPhiInCut_EB[4]    = {0.228,   0.222,   0.103,    0.0816};
  Float_t HoECut_EB[4]       = {0.356,   0.298,   0.253,    0.0414};
  Float_t diffEpCut_EB[4]    = {0.299,   0.241,   0.134,    0.0129 };
  Int_t ConvVeto_EB[4]       = {1,       1,       1,        1};
  Int_t missHitsCut_EB[4]    = {2,       1,       1,        1};

  Float_t sIeIeCut_EE[4]     = {0.037,   0.0314,  0.0298,  0.0292};
  Float_t dEtaInCut_EE[4]    = {0.00895, 0.00868, 0.00609, 0.00605};
  Float_t dPhiInCut_EE[4]    = {0.213,   0.213,   0.045,   0.0394 };
  Float_t HoECut_EE[4]       = {0.211,   0.101,   0.0878,  0.0641};
  Float_t diffEpCut_EE[4]    = {0.15,    0.14,    0.13,    0.0129 };
  Int_t ConvVeto_EE[4]       = {1,       1,       1,       1};
  Int_t missHitsCut_EE[4]    = {3,       1,       1,       1};


  if (elePt[iEle] < 10.) pass = false;
  if (fabs(eleSCEta[iEle]) > 2.5) pass = false;
  if (fabs(eleSCEta[iEle]) > 1.4442 && fabs(eleSCEta[iEle]) < 1.566) pass = false;

  if (fabs(eleSCEta[iEle]) < 1.4442) {
    if (eleSigmaIEtaIEta_Full5x5[iEle] > sIeIeCut_EB[iWP]) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) > dEtaInCut_EB[iWP]) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) > dPhiInCut_EB[iWP]) pass = false;
    if (eleHoverE[iEle] > HoECut_EB[iWP]) pass = false;
    if (eleEoverPInv[iEle] > diffEpCut_EB[iWP]) pass = false;
    if (eleConvVeto[iEle] != ConvVeto_EB[iWP]) pass = false;
    if (eleMissHits[iEle] > missHitsCut_EB[iWP]) pass = false;
  } else {
    if (eleSigmaIEtaIEta_Full5x5[iEle] > sIeIeCut_EE[iWP]) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) > dEtaInCut_EE[iWP]) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) > dPhiInCut_EE[iWP]) pass = false;
    if (eleHoverE[iEle] > HoECut_EE[iWP]) pass = false;
    if (eleEoverPInv[iEle] > diffEpCut_EE[iWP]) pass = false;
    if (eleConvVeto[iEle] != ConvVeto_EE[iWP]) pass = false;                                                                                                   
    if (eleMissHits[iEle] > missHitsCut_EE[iWP]) pass = false;
  }

  if (elePFMiniIso[iEle] > 0.1) pass = false;
  //cout << "pass or not of electron ID: " << pass << endl;

  return pass;
}

bool pass_HEEPID (int iEle) {

  bool pass = true;

  if (elePt[iEle] < 10.) pass = false;
  if (fabs(eleSCEta[iEle]) > 2.5) pass = false;
  if (fabs(eleSCEta[iEle]) > 1.4442 && fabs(eleSCEta[iEle]) < 1.566) pass = false;

  float E25OverE55 = eleE2x5Full5x5[iEle]/eleE5x5Full5x5[iEle];
  float E15OverE55 = eleE1x5Full5x5[iEle]/eleE5x5Full5x5[iEle];

  if (fabs(eleSCEta[iEle]) < 1.4442) {
    if (eleEcalDrivenSeed[iEle] != 1) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) > 0.004 ) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) >  0.06) pass = false;
    if (eleHoverE[iEle] >  1/eleEn[iEle] + 0.05  ) pass = false;
    if (fabs(eleD0[iEle]) > 0.02) pass = false;
    if (eleMissHits[iEle] > 1) pass = false;
    if (E25OverE55 <= 0.94 && E15OverE55 < 0.83) pass = false;
  }
  else {
    if (eleEcalDrivenSeed[iEle] != 1) pass = false;
    if (eleSigmaIEtaIEta_Full5x5[iEle] >0.03) pass = false;
    if (fabs(eledEtaseedAtVtx[iEle]) >  0.006) pass = false;
    if (fabs(eledPhiAtVtx[iEle]) > 0.06) pass = false;
    if (eleHoverE[iEle] >  5/eleEn[iEle] + 0.05) pass = false;
    if (fabs(eleD0[iEle]) > 0.05) pass = false;
    if (eleMissHits[iEle] > 1) pass = false;
  }

  if (elePFMiniIso[iEle] > 0.1) pass = false;

  return pass;
}

bool select_eleMVA (int iEle) {

  bool pass = true;

  if ( abs(eleSCEta[iEle]) < 0.8 && eleIDMVA[iEle] < 0.972153 ) pass = false;
  else if ( (abs(eleSCEta[iEle]) > 0.8 && abs(eleSCEta[iEle]) < 1.479) && eleIDMVA[iEle] < 0.922126) pass = false;
  else if ( abs(eleSCEta[iEle]) > 1.479 && eleIDMVA[iEle] < 0.610764) pass = false;

  return pass;
}
*/

Bool_t eleMatcher(TreeReader& data, Int_t iEle) {

  if ( !data.HasMC()) return false ;
  Int_t    nMC        = data.GetInt("nMC");
  Int_t    nEle       = data.GetInt("nEle");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Int_t*   mcMomPID   = data.GetPtrInt("mcMomPID");
  Int_t*   mcGMomPID  = data.GetPtrInt("mcGMomPID");
  Float_t* mcPt       = data.GetPtrFloat("mcPt");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  Float_t* eleEta     = data.GetPtrFloat("eleEta");
  Float_t* elePhi     = data.GetPtrFloat("elePhi");
  Float_t* mcEta      = data.GetPtrFloat("mcEta");
  Float_t* mcPhi      = data.GetPtrFloat("mcPhi");

  //return false;

  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      if ( fabs( mcPID[iMC]) != 11 ) continue ;
      if ( !(mcStatusFlag[iMC]>>0&1) || !(mcStatusFlag[iMC]>>1&1) ) continue;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],eleEta[iEle],elePhi[iEle]) < 0.1 ) return true ;
    }
  return false;

}


//Bool_t Electron_SMZg(TreeReader &data, int i) {
void Electron_SMZg(TreeReader &data, vector<int> &selectedEle) {

  Int_t nEle = data.GetInt("nEle");
  Float_t* elePt = data.GetPtrFloat("eleCalibPt");
  Float_t* eleEta = data.GetPtrFloat("eleEta");
  Float_t* eleSCEta = data.GetPtrFloat("eleSCEta");
  Float_t* eleIDMVAHZZ = data.GetPtrFloat("eleIDMVAHZZ");
  Float_t* eleD0 = data.GetPtrFloat("eleD0");
  Float_t* eleDz = data.GetPtrFloat("eleDz");
  Float_t* eleSIP = data.GetPtrFloat("eleSIP");
  Float_t* elePFChIso = data.GetPtrFloat("elePFChIso");
  Float_t* elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
  Float_t* elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
  Float_t rho = data.GetFloat("rho");

   selectedEle.clear();
   vector<int> passUnsorted;
   passUnsorted.clear();
   vector<float> elePtUnsorted;
   elePtUnsorted.clear();

   Float_t cut_d0[2] = {0.5, 0.5};
   Float_t cut_dz[2] = {1.0, 1.0};

   Int_t eta_bit = -1;

   for (int i = 0; i < nEle; i++) {
     if (elePt[i] < 10.) continue;
     if (fabs(eleSCEta[i]) > 2.5) continue;
     if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;

     if (fabs(eleSCEta[i]) >= 0.000 && fabs(eleSCEta[i]) < 1.000) eta_bit = 0;
     if (fabs(eleSCEta[i]) >= 1.000 && fabs(eleSCEta[i]) < 1.479) eta_bit = 1;
     if (fabs(eleSCEta[i]) >= 1.479 && fabs(eleSCEta[i]) < 2.000) eta_bit = 2;
     if (fabs(eleSCEta[i]) >= 2.000 && fabs(eleSCEta[i]) < 2.200) eta_bit = 3;
     if (fabs(eleSCEta[i]) >= 2.200 && fabs(eleSCEta[i]) < 2.300) eta_bit = 4;
     if (fabs(eleSCEta[i]) >= 2.300 && fabs(eleSCEta[i]) < 2.400) eta_bit = 5;
     if (fabs(eleSCEta[i]) >= 2.400 && fabs(eleSCEta[i]) < 5.000) eta_bit = 6;

     Int_t iEB = (fabs (eleSCEta[i]) <= 1.479) ? 0 : 1;

     //category for ID cut
     Int_t cat = -1;
     if (fabs(eleSCEta[i])  < 0.8)                              cat = 0;
     if (fabs(eleSCEta[i]) >= 0.8 && fabs(eleSCEta[i]) < 1.479) cat = 1;
     if (fabs(eleSCEta[i]) >= 1.479)                            cat = 2;

     //mva cut
     Float_t cut_mva_HZZ[3]    = {-0.870, -0.838, -0.763};
     Float_t lowcut_mva_HZZ[3] = {-0.211, -0.396, -0.215};

     //effective area
     Float_t cut_effA[7] = {0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};

     if (elePt[i] < 15.) continue;
     if ((elePt[i] <= 10.) && (eleIDMVAHZZ[i] < lowcut_mva_HZZ[cat])) continue;
     if (fabs(eleSCEta[i]) > 2.5) continue;
     if (fabs(eleD0[i]) > cut_d0[iEB]) continue;
     if (fabs(eleDz[i]) > cut_dz[iEB]) continue;
     if (eleSIP[i] >= 4.)              continue;

     if (elePt[i] >= 10. && eleIDMVAHZZ[i] < cut_mva_HZZ[cat])  continue;
     if (((elePFChIso[i] + std::max (0.0f, elePFNeuIso[i] + elePFPhoIso[i] - rho * cut_effA[eta_bit])) / elePt[i]) >= 0.35) continue;

     passUnsorted.push_back(i);
     elePtUnsorted.push_back(elePt[i]);
   }

   //sort electron pt in descending
   int siz = (int) passUnsorted.size();
   if (siz < 1) return;

   int ind[siz];
   TMath::Sort(siz, &elePtUnsorted.front(), ind);
   for (int i = 0; i < siz; ++i) {
     selectedEle.push_back(passUnsorted[ind[i]]);
   }
 }


void passEleId(TreeReader &data, vector<int> &selectedEle) {

  Int_t nEle = data.GetInt("nEle");
  Float_t* elePt = data.GetPtrFloat("eleCalibPt");
  Float_t* eleEta = data.GetPtrFloat("eleEta");
  Float_t* eleSCEta = data.GetPtrFloat("eleSCEta");
  Short_t* eleID = data.GetPtrShort("eleIDbit"); 

  selectedEle.clear();
  vector<int> passUnsorted;
  passUnsorted.clear();
  vector<float> elePtUnsorted;
  elePtUnsorted.clear();

  for (int i = 0; i < nEle; i++) {
    if (elePt[i] < 10.) continue;
    if (fabs(eleSCEta[i]) > 2.5) continue;
    if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
    if ( (eleID[i] >> 3 &1) == 1 ) {
      passUnsorted.push_back(i);
      elePtUnsorted.push_back(elePt[i]);
    }
  }

  //sort ele in pt descending 
  int siz = (int) passUnsorted.size();
  if (siz < 1) return;

  int ind[siz];
  TMath::Sort(siz, &elePtUnsorted.front(), ind);
  for (int i = 0; i < siz; ++i) {
    selectedEle.push_back(passUnsorted[ind[i]]);
  }

}


#endif
