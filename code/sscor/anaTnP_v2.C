#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include "untuplizer.h"
#include "tree_TnP.h"
#include "ElectronSelection.h"
#include "PhotonSelections.h"
//#include "cat.h"
#include "puweicalc.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStopwatch.h>


void anaTnP_v2(TString inputfile = "HiForest.root", TString outputfile = "XZg.root", bool mc = true)
{

  TStopwatch *countsc = new TStopwatch();
  countsc->Start();

  outputfile = "minitrees_Zee/";

  if (!mc)  outputfile = outputfile + "SingleEle";
  //else outputfile = outputfile + "data";

  if (inputfile.Contains("_Run2016B") ) outputfile += "_Run2016B";
  if (inputfile.Contains("_Run2016C") ) outputfile += "_Run2016C";
  if (inputfile.Contains("_Run2016D") ) outputfile += "_Run2016D";
  if (inputfile.Contains("_Run2016E") ) outputfile += "_Run2016E";
  if (inputfile.Contains("_Run2016F_ReMiniAOD1") ) outputfile += "_Run2016F_1";
  if (inputfile.Contains("_Run2016F_ReMiniAOD2") ) outputfile += "_Run2016F_2";
  if (inputfile.Contains("_Run2016G") ) outputfile += "_Run2016G";
  if (inputfile.Contains("_Run2016H_ReMiniAOD_03Feb_v2") ) outputfile += "_Run2016H_v2";
  if (inputfile.Contains("_Run2016H_ReMiniAOD_03Feb_v3") ) outputfile += "_Run2016H_v3";
  if (inputfile.Contains("_Zg_aMCatNLO/") ) outputfile += "_Zg_aMCatNLO";
  if (inputfile.Contains("_Zg_pt130/") ) outputfile += "_Zg_pt130";
  if (inputfile.Contains("DYJetsToLL_m50_MG") ) outputfile += "DYJets_MG";
  if (inputfile.Contains("DYJetsToLL_m50_aMCatNLO") ) outputfile += "DYJets_amcatnlo";
  if (inputfile.Contains("WJetsToLNu") ) outputfile += "WJetsToLNu";
  if (inputfile.Contains("TT") ) outputfile += "TTbar";


  outputfile += "_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr.root";
  //outputfile += "_TnP_Zee_BDT_Upto6000_Corr_sieie_sieip_etawidth_r9_s4_phiwidth.root";


  TreeReader data(inputfile, "ggNtuplizer/EventTree");

  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *passingIdTree = new TTree("passingIdTree", "passing Id tree");

  inittree(passingIdTree);

  Long64_t nTotal(0) ;
  Long64_t nPassHLT(0) ;
  Long64_t nPassTag(0) ;
  Long64_t nFill(0) ;
  Long64_t ntot(0);

  PUWeightCalculator puCalcGJ;
  PUWeightCalculator puCalcGJ_69nb;
  PUWeightCalculator puCalcGJ_69p2nb;
  PUWeightCalculator puCalcGJ_65nb;
  PUWeightCalculator puCalcGJ_63nb;

  if ( data.HasMC() ) {
    //puCalcGJ.Init("external/PU_histo_13TeV_GoldenJSON_summer16_69000nb.root");
    puCalcGJ_69nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_72400nb.root");
    puCalcGJ_69p2nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_69200nb.root");
    puCalcGJ_65nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_66000nb.root");
    puCalcGJ_63nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_63000nb.root");
  }


  //for scale factors
  TFile *feleRecoSF = new TFile("external/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root", "READ");
  TFile *feleSF = new TFile("external/2016LegacyReReco_ElectronTight_Fall17V2.root", "READ");
  TFile *fphoSF = new TFile("external/phoSelSF_Zg.root", "READ");

  TH2F *heleRecoSF = (TH2F*) feleRecoSF->Get("EGamma_SF2D");
  TH2F *heleSF = (TH2F*) feleSF->Get("EGamma_SF2D");
  TH2F *hphoSF = (TH2F*) fphoSF->Get("EGamma_SF2D");

  //showershape correction
  TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19_sieieCorr.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19_sieie_sieip_Corr.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19_sieie_sieip_etawidth_Corr.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19_sieie_sieip_etawidth_r9_Corr.root", "read");
  //TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000_May19_sieie_sieip_etawidth_r9_s4_Corr.root", "read");

  TGraph *tgr[14];
  tgr[0] = (TGraph*) fss_rw->Get("transfEtaWidthEB");
  tgr[1] = (TGraph*) fss_rw->Get("transfS4EB");
  tgr[2] = (TGraph*) fss_rw->Get("transffull5x5R9EB");
  tgr[3] = (TGraph*) fss_rw->Get("transffull5x5sieieEB");
  tgr[4] = (TGraph*) fss_rw->Get("transffull5x5sieipEB");
  tgr[5] = (TGraph*) fss_rw->Get("transfrhoEB");
  tgr[6] = (TGraph*) fss_rw->Get("transfPhiWidthEB");

  //tgr[4] = (TGraph*) fss_rw->Get("transfEtaWidthEE");
  //tgr[5] = (TGraph*) fss_rw->Get("transfS4EE");
  //tgr[6] = (TGraph*) fss_rw->Get("transffull5x5R9EE");
  //tgr[7] = (TGraph*) fss_rw->Get("transffull5x5sieieEE");
  tgr[7] = (TGraph*) fss_rw->Get("transfEtaWidthEE");
  tgr[8] = (TGraph*) fss_rw->Get("transfS4EE");
  tgr[9] = (TGraph*) fss_rw->Get("transffull5x5R9EE");
  tgr[10] = (TGraph*) fss_rw->Get("transffull5x5sieieEE");
  tgr[11] = (TGraph*) fss_rw->Get("transffull5x5sieipEE");
  tgr[12] = (TGraph*) fss_rw->Get("transfrhoEE");
  tgr[13] = (TGraph*) fss_rw->Get("transfPhiWidthEE");
  
  Long64_t TotalEntries = data.GetEntriesFast();  
  for ( Long64_t ev = 0; ev < TotalEntries; ev++) {   
    //for ( Long64_t ev = 0; ev < 10000000; ev++) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);
    readggtree(data);

    nTotal++ ;

    if ( nEle < 1 ) continue ;

    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    isPVGood_ = isPVGood;
    run_ = run;
    event_ = event;
    lumi_ = lumis;
    rho_ = rho;

    //pfmet
    pfMet = pfMET;
    pfMetPhi = pfMETPhi;

    //PU weight for MC
    if (data.HasMC()) {
      float* puTrue = data.GetPtrFloat("puTrue");
      //puweigj = (float) puCalcGJ.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69nb = (float) puCalcGJ_69nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69p2nb = (float) puCalcGJ_69p2nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_65nb = (float) puCalcGJ_65nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_63nb = (float) puCalcGJ_63nb.GetWeight(run, puTrue[1]); // in-time PU

      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    //trigger
    trig_Ele27_Eta2p1_WPTight = hlt>>1&1;
    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele23_Ele12 = hlt>>5&1;
    trig_DoubleEle33 = hlt>>10&1;

    if (!trig_Ele27_WPTight) continue;
    
    vector<int> acc_tag_ele, acc_probe_ele;
    int nZ = 0;

    for (int i=0; i<nEle; i++) {
      if (elePt[i] < 25.) continue;
      if (fabs(eleSCEta[i]) > 2.5) continue;
      //if (fabs(eleSCEta[i]) > 1.4442 && fabs(eleSCEta[i]) < 1.566) continue;
      if (mc) {
        bool match = eleMatcher(data, i);
        if (!match) continue;
      }

      if ( (eleID[i] >>3&1)==1 && (eleFiredSingleTrgs[i] >>12&1)==1) acc_tag_ele.push_back(i);
      acc_probe_ele.push_back(i);
    }

    if (acc_tag_ele.size() < 1) continue;

    nTag = 0;

    for (vector<int>::iterator itag = acc_tag_ele.begin(); itag != acc_tag_ele.end(); itag++) {
      if (fabs(eleSCEta[*itag]) > 1.4442 && fabs(eleSCEta[*itag]) < 1.566) continue;
      if (fabs(eleSCEta[*itag]) > 2.1) continue;
      if ( elePt[*itag] < 30.) continue;

      //initial value;
      Tag_Pt = -99.;
      Tag_Eta = -99.;
      Tag_Phi = -99.;
      Tag_SCEta = -99.;
      TagIndex = -99;
      Tag_sieie = -99.;
      Tag_chIso = -99.;
      Tag_neuIso = -99.;
      Tag_phoIso = -99.;
      Tag_hoe = -99.;
      Tag_eop = -99.;
      Tag_ch = -99.;
      Tag_dEtaSeed = -99.;
      Tag_dEtaVtx = -99.;
      Tag_dPhiVtx = -99.;

      nProbe = 0;
      nPassId = 0;
      nFailId = 0;
      nPassIdIso = 0;
      nFailIdIso = 0;
      ProbeIndex = -99;
      Probe_Pt = -99.;
      Probe_Eta = -99.;
      Probe_Phi = -99.;
      Probe_SCEta = -99.;
      ProbeIndex = -99;
      Probe_sieie = -99.;
      Probe_chIso = -99.;
      Probe_neuIso = -99.;
      Probe_phoIso = -99.;
      Probe_hoe = -99.;
      Probe_eop = -99.;
      Probe_ch = -99.;
      Probe_dEtaSeed = -99.;
      Probe_dEtaVtx = -99.;
      Probe_dPhiVtx = -99.;
      Probe_miniIso = -99.;
      LooseID = 0;
      MediumID = 0;
      TightID = 0;
      LooseIDOnly = 0;
      Probe_PhoEleVeto = 0;
      Probe_PhoChIso = -99.;
      passPhoID_Zg = 0;
      Probe_ssmva = -99.;
      Probe_ssmva_rw = -99.;
      pho_presel = 0;

      Zm = -99.;
      Zeta = -99.;
      Zpt = -99.;

      TLorentzVector vTag ;
      vTag.SetPtEtaPhiM(elePt[*itag], eleEta[*itag], elePhi[*itag], 0.511*0.001);

      vector<float> Zmass, ZEta, ZPt;
      vector<float> ProbePt;
      vector<int> index;
      nProbe = 0;

      //for (vector<int>::iterator iele = acc_probe_ele.begin(); iele != acc_probe_ele.end(); iele++) {
      //if ( *itag == *iele) continue;
      //TLorentzVector vProbe ,vZ ;

      //vProbe.SetPtEtaPhiM(elePt[*iele], eleEta[*iele], elePhi[*iele],0.511*0.001);
      //vZ = vTag + vProbe;

      if (nPho < 1) continue;
      for (int ipho = 0; ipho < nPho; ipho++) {
	if (phoEt[ipho] < 15.) continue;
	if (fabs(phoSCEta[ipho]) > 2.5) continue;
	if (fabs(phoSCEta[ipho]) > 1.4442 && fabs(phoSCEta[ipho]) < 1.566) continue;

	//not match to tag ele
	if ( fabs(phoSCEta[ipho] - eleSCEta[*itag]) <= 0.05 || fabs(deltaPhi(phoSCPhi[ipho], eleSCPhi[*itag]) <= 0.05) ) continue;

	//match to gen ele
	//if ( data.HasMC()) 
	  //if (phoMatcher(data, ipho) == 1) continue; //for WJets

	TLorentzVector vProbe ,vZ ;
	vProbe.SetPtEtaPhiM(phoEt[ipho] , phoEta[ipho] , phoPhi[ipho] , 0. );
	vZ = vTag + vProbe; 

	if (mc) {
	  if ( vZ.M() < 40 || vZ.M() > 140 ) continue ;  //for MC template
	} else {
	  if ( vZ.M() < 60 || vZ.M() > 120 ) continue ; //for data
	}
	  
	
	Zmass.push_back(vZ.M());
	ZEta.push_back(vZ.Eta());
	ZPt.push_back(vZ.Pt());
	//ProbePt.push_back(elePt[*iele]);
	//index.push_back(*iele);
	ProbePt.push_back(phoEt[ipho]);
	index.push_back(ipho);
      }
    
      nZ = Zmass.size();
      if (nZ < 1) continue;
      float minZ = 999.;
      int iProbe = -1;

      //choose Z having closest to truth mass
      /*for (int i = 0; i < nZ; i++) {
	//cout << "Zmass: " << Zmass[i];
	float deltaM = fabs(91.19-Zmass[i]);
	if (minZ > deltaM) {
	  minZ = deltaM;
	  iProbe = index[i];
	  Zm = Zmass[i];
	  Zeta = ZEta[i];
	  Zpt = ZPt[i];
	}
	}*/

      //choose highest probe pt
      float maxpt = 0.;
      for (int i = 0; i < nZ; i++) {
	if (ProbePt[i] > maxpt) {
	  maxpt = ProbePt[i];
	  Zm = Zmass[i];
          Zeta = ZEta[i];
          Zpt = ZPt[i];
	  iProbe = index[i];
	}
      }

    
      if (iProbe < 0) continue;
      //if ( *itag == iProbe) continue;

      bool passId = false;
      bool passIdIso = false;

      nTag++;    
      Tag_Pt = elePt[*itag];
      Tag_Eta = eleEta[*itag];
      Tag_SCEta= eleSCEta[*itag];
      Tag_Phi = elePhi[*itag];
      Tag_sieie = eleSigmaIEtaIEta_Full5x5[*itag];
      Tag_chIso = elePFChIso[*itag];
      Tag_neuIso = elePFNeuIso[*itag];
      Tag_phoIso = elePFPhoIso[*itag];
      Tag_hoe = eleHoverE[*itag];
      Tag_eop = eleEoverP[*itag];
      //Tag_dEtaSeed = eledEtaseedAtVtx[*itag];
      Tag_dEtaVtx = eledEtaAtVtx[*itag];
      Tag_dPhiVtx = eledPhiAtVtx[*itag];
      Tag_ch = eleCharge[*itag];
      //Tag_IDMVA = eleIDMVA[*itag];
      ele1_MatchTrig27 = eleFiredSingleTrgs[*itag] >>12&1;
      mW_Tag_Met = sqrt(2*pfMet*elePt[*itag]*(1-cos(deltaPhi(pfMetPhi,elePhi[*itag]))));

      //scale factor
      float eleRecoSF = 1.;
      float eleSF = 1.;
      if ( data.HasMC() ) {
        //reco sf
	float recopt = elePt[*itag];
	int NBinsYRec = heleRecoSF->GetNbinsY();
	float maxPt = heleRecoSF->GetYaxis()->GetBinCenter(NBinsYRec);
	if ( elePt[*itag] > maxPt) recopt = maxPt;
	if (elePt[*itag] < 25.) recopt = 25.;

	int binyrec = heleRecoSF->GetYaxis()->FindBin(recopt);
        int binxrec = heleRecoSF->GetXaxis()->FindBin(eleSCEta[*itag]);
        eleRecoSF = heleRecoSF->GetBinContent(binxrec, binyrec);

        //ID SFs
	float selpt = elePt[*itag];
	int NBinsY = heleSF->GetNbinsY();
	maxPt = heleSF->GetYaxis()->GetBinCenter(NBinsY);
	if (selpt > maxPt) selpt = maxPt;
        int biny = heleSF->GetYaxis()->FindBin(selpt);
        int binx = heleSF->GetXaxis()->FindBin(eleSCEta[*itag]);
	eleSF = heleSF->GetBinContent(binx, biny);
      }
      Tag_RecoSF = eleRecoSF;
      Tag_SelSF = eleSF;


      /*
      Probe_Pt = elePt[iProbe];
      Probe_Eta = eleEta[iProbe];
      Probe_SCEta = eleSCEta[iProbe];
      Probe_Phi = elePhi[iProbe];
      Probe_sieie = eleSigmaIEtaIEta_Full5x5[iProbe];
      Probe_chIso = elePFChIso[iProbe];
      Probe_neuIso = elePFNeuIso[iProbe];
      Probe_phoIso = elePFPhoIso[iProbe];
      Probe_hoe = eleHoverE[iProbe];
      Probe_eop = eleEoverP[iProbe];
      Probe_dEtaSeed = eledEtaseedAtVtx[iProbe];
      Probe_dEtaVtx = eledEtaAtVtx[iProbe];
      Probe_dPhiVtx = eledPhiAtVtx[iProbe];
      Probe_ch = eleCharge[iProbe];
      Probe_miniIso = elePFMiniIso[iProbe];

      LooseID = eleID[iProbe] >>1&1;
      MediumID = eleID[iProbe] >>2&1;
      TightID = eleID[iProbe] >>3&1;
      if ( pass_cutbased_80X(iProbe,0) ) LooseIDOnly = 1;
      else LooseIDOnly = 0;
      */
      //cout << "filling probe" << endl;
      Probe_Pt = phoEt[iProbe];
      Probe_Eta = phoEta[iProbe];
      Probe_SCEta = phoSCEta[iProbe];
      Probe_Phi = phoPhi[iProbe];
      Probe_PhoEleVeto = phoEleVeto[iProbe];
      //Probe_PhoChIso = phoPFChIso[iProbe];
      Probe_ssmva = PhotonSSMVA(data, iProbe);
      Probe_ssmva_rw = PhotonSSMVA_rw(data, iProbe, tgr);

      Probe_hoe = phoHoverE[iProbe];
      Probe_ChWorstIso = phoPFChWorstIso[iProbe];
      Probe_SCEtaWidth = phoSCEtaWidth[iProbe];
      Probe_s4Full5x5 = phoE2x2Full5x5[iProbe]/phoE5x5Full5x5[iProbe];
      Probe_R9 = phoR9[iProbe];
      Probe_sieie = phoSigmaIEtaIEtaFull5x5[iProbe];
      Probe_sieip = phoSigmaIEtaIPhiFull5x5[iProbe];
      Probe_SCRawE = phoSCRawE[iProbe];
      Probe_scphiwidth = phoSCPhiWidth[iProbe];


      if (data.HasMC()) {
        if(TMath::Abs(phoSCEta[iProbe])<1.5) {
          //Probe_SCEtaWidth_rw    = tgr[3]->Eval(Probe_SCEtaWidth);
	  //float etaw    = tgr[3]->Eval(Probe_SCEtaWidth);
          //Probe_SCEtaWidth_rw    = tgr[0]->Eval(etaw);

          Probe_SCEtaWidth_rw    = tgr[0]->Eval(Probe_SCEtaWidth);
          Probe_s4Full5x5_rw     = tgr[1]->Eval(Probe_s4Full5x5);
          Probe_R9_rw            = tgr[2]->Eval(Probe_R9);
          Probe_sieie_rw         = tgr[3]->Eval(Probe_sieie);
	  Probe_sieip_rw         = tgr[4]->Eval(Probe_sieip);
	  rho_rw                 = tgr[5]->Eval(rho);
	  Probe_SCPhiWidth_rw    = tgr[6]->Eval(Probe_scphiwidth);
        } else {
          //Probe_SCEtaWidth_rw    = tgr[10]->Eval(Probe_SCEtaWidth);
	  //float etaw    = tgr[10]->Eval(Probe_SCEtaWidth);
          //Probe_SCEtaWidth_rw    = tgr[7]->Eval(etaw);

          Probe_SCEtaWidth_rw    = tgr[7]->Eval(Probe_SCEtaWidth);
          Probe_s4Full5x5_rw     = tgr[8]->Eval(Probe_s4Full5x5);
          Probe_R9_rw            = tgr[9]->Eval(Probe_R9);
          Probe_sieie_rw         = tgr[10]->Eval(Probe_sieie);
	  Probe_sieip_rw         = tgr[11]->Eval(Probe_sieip);
          rho_rw                 = tgr[12]->Eval(rho);
          Probe_SCPhiWidth_rw    = tgr[13]->Eval(Probe_scphiwidth);

        }
      }


      //rho-correction charged isolation 
      float chEA[]  = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
      float neuEA[] = {0.0597, 0.0807, 0.0629, 0.0197, 0.0184, 0.0284, 0.0591};
      float phoEA[] = {0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998};
      int bin;
      if (phoSCEta[iProbe] < 1.0) bin = 0;
      else if (phoSCEta[iProbe] >= 1.0 && phoSCEta[iProbe] < 1.479) bin = 1;
      else if (phoSCEta[iProbe] >= 1.479 && phoSCEta[iProbe] < 2.0) bin = 2;
      else if (phoSCEta[iProbe] >= 2.0 && phoSCEta[iProbe] < 2.2) bin = 3;
      else if (phoSCEta[iProbe] >= 2.2 && phoSCEta[iProbe] < 2.3) bin = 4;
      else if (phoSCEta[iProbe] >= 2.3 && phoSCEta[iProbe] < 2.4) bin = 5;
      else bin = 6; // phoSCEta[iProbe] >= 2.4

      double corr = phoPFChIso[iProbe] - chEA[bin] * rho;
      double corr_neu = phoPFNeuIso[iProbe] - neuEA[bin] * rho;
      double corr_pho = phoPFPhoIso[iProbe] - phoEA[bin] * rho;

      double corrchIso = TMath::Max(corr, 0.);
      double corrNeuIso = TMath::Max(corr_neu, 0.);
      double corrPhoIso = TMath::Max(corr_pho, 0.);

      Probe_PhoChIso = corrchIso;
      Probe_neuIso = corrNeuIso;
      Probe_phoIso = corrPhoIso;

      //loose cut-based ID
      if ( fabs(phoSCEta[iProbe]) < 1.5) {
	if (Probe_hoe < 0.0597 && Probe_sieie < 0.01031 && Probe_PhoChIso < 1.295) {
	  if (Probe_neuIso < 10.910 + 0.0148*Probe_Pt + 0.000017*Probe_Pt*Probe_Pt) {
	    if ( Probe_phoIso < 3.630 + 0.0047*Probe_Pt) 
	      LooseID = 1;
	  }
	}
      }
      else if ( fabs(phoSCEta[iProbe]) > 1.5) {
	if (Probe_hoe < 0.0481 && Probe_sieie < 0.03013 && Probe_PhoChIso < 1.011) {
          if (Probe_neuIso < 5.931 + 0.0163*Probe_Pt + 0.000014*Probe_Pt*Probe_Pt) {
            if ( Probe_phoIso < 6.641 + 0.0034*Probe_Pt)
              LooseID = 1;
          }
        }
      }

      pho_presel = preselect_pho_Zg(iProbe);
      if ( preselect_pho_Zg(iProbe) == 1 )
	if ( (fabs(phoSCEta[iProbe]) < 1.5 && Probe_PhoChIso<2.) || (fabs(phoSCEta[iProbe]) > 1.5 && Probe_PhoChIso <1.5) )
	  passPhoID_Zg = 1;
     

      //photon scale factor
      float phoSF = 1., phoCSEV_SF = 1.;
      if ( mc ) {
        int binx, biny;
        biny = hphoSF->GetYaxis()->FindBin(phoEt[iProbe]);
	binx = hphoSF->GetXaxis()->FindBin(phoSCEta[iProbe]);
        phoSF = hphoSF->GetBinContent(binx, biny);
        //phoCSEV_SF = hphoCSEV->GetBinContent(binx, 1);
        int NBinsY = hphoSF->GetNbinsY();
        float maxEt = hphoSF->GetYaxis()->GetBinUpEdge(NBinsY);
        if ( phoCalibEt[iProbe] > maxEt) phoSF = hphoSF->GetBinContent(binx, NBinsY);

      } 
      Probe_phoSF = phoSF;

      passingIdTree->Fill();

    }

  }// end of loop
  cout << "total events = " << nTotal << endl ;
  cout << "pass HLT     = " << nPassHLT << endl ;
  cout << "pas tags     = " << nPassTag << endl ;
  cout << "filled       = " << nFill << endl ;
  cout << "endloop " << endl ;
  
  fo->Write();
  fo->Close();
  delete fo;

  countsc->Stop();
  countsc->Print();


}
