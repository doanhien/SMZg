#ifndef MuonSelection_h__
#define MuonSelection_h__

//#include "xAna.h"

void Muon_SMZg(TreeReader &data, vector<int> &selectedMu) {

  Int_t nMu = data.GetInt("nMu");
  Float_t* muPt = data.GetPtrFloat("muPt");
  Float_t* muEta = data.GetPtrFloat("muEta");
  Float_t* muD0 = data.GetPtrFloat("muD0");
  Float_t* muDz = data.GetPtrFloat("muDz");
  Int_t* muBestTrkType = data.GetPtrInt("muBestTrkType");
  Int_t* muType = data.GetPtrInt("muType");
  Float_t* muBestTrkPtError = data.GetPtrFloat("muBestTrkPtError");
  Float_t* muBestTrkPt = data.GetPtrFloat("muBestTrkPt");
  Int_t* muStations = data.GetPtrInt("muStations");
  Int_t* muPixelHits = data.GetPtrInt("muPixelHits");
  Int_t* muTrkLayers = data.GetPtrInt("muTrkLayers");
  Float_t* muSIP = data.GetPtrFloat("muSIP");
  Float_t* muPFChIso03 = data.GetPtrFloat("muPFChIso03");
  Float_t* muPFPhoIso03 = data.GetPtrFloat("muPFPhoIso03");
  Float_t* muPFNeuIso03 = data.GetPtrFloat("muPFNeuIso03");
  Float_t* muPFPUIso03 = data.GetPtrFloat("muPFPUIso03");

  selectedMu.clear();
  vector<int> passUnsorted;
  passUnsorted.clear();
  vector<float> muPtUnsorted;
  muPtUnsorted.clear();

  for (int i = 0; i < nMu; i++) {
    if (muPt[i] < 10.) continue;
    if (fabs(muEta[i]) > 2.4) continue;

    if (muSIP[i]       > 4.)                                                           continue;
    if (fabs(muD0[i])  > 0.5)                                                          continue;
    if (fabs(muDz[i])  > 1.)                                                           continue;
    if (muBestTrkType[i] == 2)                                                         continue;
    if ((muPt[i]<200) && ((muType[i] >> 5 & 1) == 0) )                                 continue;
    if (((muType[i]>>1&1)==0) && (((muType[i] >> 2 & 1) == 0) && (muStations[i] >0)))  continue;
    if (muPt[i] > 200.)
      {
        if ((muType[i] >> 5 & 1)== 0)
          {
            if ((muType[i] >> 2 & 1) == 0)                continue;
            if (muBestTrkPtError[i]/muBestTrkPt[i]>=0.3)  continue;
            if (muStations[i]   < 2)    continue;
            if (fabs(muD0[i])  >= 0.2)  continue;
            if (fabs(muDz[i])  >= 0.5)  continue;
            if (muPixelHits[i] <= 0)    continue;
            if (muTrkLayers[i] <= 5)    continue;
          }
      }
    if ((muPFChIso03[i] + TMath::Max(0., muPFPhoIso03[i]+muPFNeuIso03[i]-0.5*muPFPUIso03[i])) / muPt[i] > 0.35)  continue;

    passUnsorted.push_back(i);
    muPtUnsorted.push_back(muPt[i]);
  }

  //sort pt in descending
  int siz = (int) passUnsorted.size();
  if (siz < 1) return;

  int ind[siz];
  TMath::Sort(siz, &muPtUnsorted.front(), ind);
  for (int i = 0; i < siz; ++i) {
    selectedMu.push_back(passUnsorted[ind[i]]);
  }

}

void passMuonId (TreeReader &data, vector<int> &selectedMu) {

  Int_t nMu = data.GetInt("nMu");
  Float_t* muPt = data.GetPtrFloat("muPt");
  Float_t* muEta = data.GetPtrFloat("muEta");
  Short_t* muIDbit = data.GetPtrShort("muIDbit");
  Float_t* muPFChIso = data.GetPtrFloat("muPFChIso");
  Float_t* muPFPhoIso = data.GetPtrFloat("muPFPhoIso");
  Float_t* muPFNeuIso = data.GetPtrFloat("muPFNeuIso");
  Float_t* muPFPUIso = data.GetPtrFloat("muPFPUIso");

  selectedMu.clear();
  vector<int> passUnsorted;
  passUnsorted.clear();
  vector<float> muPtUnsorted;
  muPtUnsorted.clear();

  for (int i = 0; i < nMu; i++) {
    if (muPt[i] < 10.) continue;
    if (fabs(muEta[i]) > 2.4) continue;
    if ( ((muIDbit[i] >> 2) &1) != 1) continue;
    if ( (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] > 0.15) continue;
    passUnsorted.push_back(i);
    muPtUnsorted.push_back(muPt[i]);
  }

  //sort pt in descending
  int siz = (int) passUnsorted.size();
  if (siz < 1) return;

  int ind[siz];
  TMath::Sort(siz, &muPtUnsorted.front(), ind);
  for (int i = 0; i < siz; ++i) {
    selectedMu.push_back(passUnsorted[ind[i]]);
  }

}


Bool_t muonMatcher(TreeReader& data, Int_t iMu) {

  if ( !data.HasMC()) return false ;
  Int_t    nMC        = data.GetInt("nMC");
  Int_t*   mcPID      = data.GetPtrInt("mcPID");
  Int_t*   mcMomPID   = data.GetPtrInt("mcMomPID");
  Int_t*   mcGMomPID  = data.GetPtrInt("mcGMomPID");
  Float_t* mcPt       = data.GetPtrFloat("mcPt");
  Float_t* mcEta       = data.GetPtrFloat("mcEta");
  Float_t* mcPhi       = data.GetPtrFloat("mcPhi");
  UShort_t* mcStatusFlag = (UShort_t*) data.GetPtrShort("mcStatusFlag");
  Float_t* muPt = data.GetPtrFloat("muPt");
  Float_t* muEta = data.GetPtrFloat("muEta");
  Float_t* muPhi = data.GetPtrFloat("muPhi");


  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      if ( fabs( mcPID[iMC]) != 13 ) continue ;
      //if (       mcMomPID[iMC]  != 23 && mcMomPID[iMC]  != 24 ) continue ;
      if ( ((mcStatusFlag[iMC]>>0)&1)!=1 || ((mcStatusFlag[iMC]>>1)&1)!=1 ) continue;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],muEta[iMu],muPhi[iMu]) < 0.1 ) return true ;
    }
  return false;
}


#endif
