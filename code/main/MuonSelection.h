#ifndef MuonSelection_h__
#define MuonSelection_h__

//#include "xAna.h"

bool LooseMuon (int imu) {
  if ( ((muType[imu]>>5) & 1)==1
       && ( (muType[imu]>>1 & 1)==1 || (muType[imu]>>2 & 1)==1) 
       && muPFMiniIso[imu] < 0.1
       ) return true;
  else return false;
}

bool HighPtMuon (int imu) {
  if ( (muType[imu]>>1 & 1)==1 && (muType[imu]>>2 & 1)==1
       && muMuonHits[imu]>0
       && muPixelHits[imu]>0 
       && muStations[imu]>1
       && fabs(muInnerD0[imu])<0.2
       && fabs(muInnerDz[imu])<0.5
       && muTrkLayers[imu]>5
       && muBestTrkPtError[imu]/muBestTrkPt[imu] < 0.3
       && muPFMiniIso[imu] < 0.1
       ) return true;  
  else return false;


}

bool goodMuon (int imu) {

  if(muChi2NDF[imu]<10
     && muMuonHits[imu]>0
     && muStations[imu]>1
     && muTrkLayers[imu]>5
     && muPixelHits[imu]>0
     && fabs(muInnerD0[imu])<0.2
     && fabs(muInnerDz[imu])<0.5
      ) return true;
  else return false;
}

/*
Bool_t Muon_SMZg(int i) {

  if (muPt[i] < 5. )         return false;  
  if (fabs(muEta[i]) > 2.4)  return false;

  if (muSIP[i]       > 4.)                                                           return false;
  if (fabs(muD0[i])  > 0.5)                                                          return false;
  if (fabs(muDz[i])  > 1.)                                                           return false;
  if (muBestTrkType[i] == 2)                                                         return false;
  if ((muPt[i]<200) && ((muType[i] >> 5 & 1) == 0) )                                 return false;
  if (((muType[i]>>1&1)==0) && (((muType[i] >> 2 & 1) == 0) && (muStations[i] >0)))  return false;
  if (muPt[i] > 200.)
    {
      if ((muType[i] >> 5 & 1)== 0)
	{
	  if ((muType[i] >> 2 & 1) == 0)                return false;
	  if (muBestTrkPtError[i]/muBestTrkPt[i]>=0.3)  return false;
	  if (muStations[i]   < 2)    return false;
	  if (fabs(muD0[i])  >= 0.2)  return false;
	  if (fabs(muDz[i])  >= 0.5)  return false;
	  if (muPixelHits[i] <= 0)    return false;
	  if (muTrkLayers[i] <= 5)    return false;
	}
    }
  if ((muPFChIso03[i] + TMath::Max(0., muPFPhoIso03[i]+muPFNeuIso03[i]-0.5*muPFPUIso03[i])) / muPt[i] > 0.35)  return false;

  return true;


  }*/


void Muon_SMZg (vector<int> &selectedMu) {
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


void passMuonId (vector<int> &selectedMu) {
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
    //if ( (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] > 0.25) continue;
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

void NoMuonId (vector<int> &selectedMu) {
  selectedMu.clear();
  vector<int> passUnsorted;
  passUnsorted.clear();
  vector<float> muPtUnsorted;
  muPtUnsorted.clear();

  for (int i = 0; i < nMu; i++) {
    if (muPt[i] < 10.) continue;
    if (fabs(muEta[i]) > 2.4) continue;
    //if ( ((muIDbit[i] >> 2) &1) != 1) continue;
    //if ( (muPFChIso[i] + TMath::Max(0., muPFPhoIso[i] + muPFNeuIso[i] -0.5 * muPFPUIso[i]))/muPt[i] > 0.15) continue;
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


  for (int iMC = 0; iMC < nMC; ++iMC)
    {
      //if (mcPt[iMC]< 10.) continue ;
      //if (fabs(mcEta[iMC]) > 2.5) continue;
      if ( fabs( mcPID[iMC]) != 13 ) continue ;
      if ( ((mcStatusFlag[iMC]>>0)&1)!=1 || ((mcStatusFlag[iMC]>>1)&1)!=1 ) continue;
      //if ( ((mcStatusFlag[iMC]>>0)&1) == 0 ) continue;
      if ( deltaR(mcEta[iMC],mcPhi[iMC],muEta[iMu],muPhi[iMu]) < 0.1 ) return true ;
    }
  return false;
}

#endif
