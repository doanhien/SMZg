#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TF1.h>
#include <TMath.h>

#include <TRandom3.h>

#include "untuplizer.h"
#include "tree.h"
#include "PhotonSelections.h"
#include "MuonSelection.h"
#include "ElectronSelection.h"
#include "puweicalc.h"
#include "getSFs.h"

#include "roccor.Run2.v3/RoccoR.cc"

void xAna(TString inputfile = "HiForest.root", TString outputfile = "Zg.root", bool mc = true, int iWP = 0) {


  outputfile = "minitrees/";
  if (inputfile.Contains("SingleEle") ) outputfile += "SingleEle";
  else if (inputfile.Contains("DoubleEG") ) outputfile += "DoubleEG";
  else if (inputfile.Contains("SingleMuon") ) outputfile += "SingleMuon";
  else if (inputfile.Contains("DoubleMu") ) outputfile += "DoubleMu";
  else outputfile += "";


  if (inputfile.Contains("_Run2016B_Legacy") ) outputfile += "_Run2016B_Legacy";
  if (inputfile.Contains("_Run2016C_Legacy") ) outputfile += "_Run2016C_Legacy";
  if (inputfile.Contains("_Run2016D_Legacy") ) outputfile += "_Run2016D_Legacy";
  if (inputfile.Contains("_Run2016E_Legacy") ) outputfile += "_Run2016E_Legacy";
  if (inputfile.Contains("_Run2016F_Legacy") ) outputfile += "_Run2016F_Legacy";
  if (inputfile.Contains("_Run2016G_Legacy") ) outputfile += "_Run2016G_Legacy";
  if (inputfile.Contains("_Run2016H_Legacy") ) outputfile += "_Run2016H_Legacy";

  if (inputfile.Contains("_Zg_aMCatNLO") ) outputfile += "Zg_aMCatNLO";
  if (inputfile.Contains("_Zg_pt130") ) outputfile += "Zg_pt130";
  if (inputfile.Contains("_DYJetsToLL_m50_aMCatNLO") ) outputfile += "ZJets_aMCatNLO";
  if (inputfile.Contains("_DYJetsToLL_m50_MG") ) outputfile += "ZJets_MG";
  if (inputfile.Contains("_TTTo2L2Nu_powheg") ) outputfile += "TT_Powheg";
  if (inputfile.Contains("_ZZTo4L") ) outputfile += "ZZTo4L";
  if (inputfile.Contains("_WZTo3LNu") ) outputfile += "WZTo3LNu";
  if (inputfile.Contains("_WWTo2L2Nu") ) outputfile += "WWTo2L2Nu";
  if (inputfile.Contains("_WWToLNuQQ") ) outputfile += "WWToLNuQQ";
  if (inputfile.Contains("_ZZTo2L2Nu") ) outputfile += "ZZTo2L2Nu";
  if (inputfile.Contains("_ZZTo2L2Q") ) outputfile += "ZZTo2L2Q";
  if (inputfile.Contains("_WZTo2L2Q") ) outputfile += "WZTo2L2Q";
  if (inputfile.Contains("ZGToLLG_5f") ) outputfile += "ZGToLLG_5f";
  if (inputfile.Contains("ZGToLLG_0123j") ) outputfile += "ZGToLLG_LO_Madgraph";
  if (inputfile.Contains("sherpa") ) outputfile += "ZGToLLG_5f_sherpa";


  //outputfile += "_Summer16_TMVA420_UpTp6000_5VarCorr_isPVGood_allSFs_NoGenPtCut.root";
  //outputfile += "_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  outputfile += "_NoPtCut.root";

  TreeReader data1(inputfile, "ggNtuplizer/EventTree");

  //print type of variables
  //data1.Print();


  TFile *fo = TFile::Open(outputfile, "RECREATE");

  TTree *outtree = new TTree("outtree", "output tree");
  //TTree *tZ = new TTree("tZ", "output tree");
  inittree(outtree);
  //inittree(tZ);

  TTree *outtreeGen = new TTree("outtreeGen", "output tree at gen level");
  inittreeGen(outtreeGen);

  TH1F *heemass = new TH1F("heemass", "; m^{ee} (GeV/c^{2}); Events", 30, 60, 120);
  TH1F *heemassSS = new TH1F("heemassSS", "; m^{ee} (GeV/c^{2}); Events", 30, 60, 120);
  TH1D *hntotweight = new TH1D("hntotweight", "", 2, 1, 3);

  Long64_t ev1 = 0;
  Long64_t ev2 = -1;

  if (ev2 < 0) ev2 = data1.GetEntriesFast();
  if (ev2 > data1.GetEntriesFast()) ev2 = data1.GetEntriesFast();

  RoccoR  rc("roccor.Run2.v3/RoccoR2016.txt");

  //pile up reweighting for MC
  PUWeightCalculator puCalcGJ_69nb;
  PUWeightCalculator puCalcGJ_69p2nb;
  PUWeightCalculator puCalcGJ_65nb;
  PUWeightCalculator puCalcGJ_63nb;

  //PUWeightCalculator puCalcSJ;
  if ( data1.HasMC() ) {
    puCalcGJ_69nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_72400nb.root");
    puCalcGJ_69p2nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_69200nb.root");
    puCalcGJ_65nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_66000nb.root");
    puCalcGJ_63nb.Init("external/PU_histo_13TeV_2016_GoldenJSON_63000nb.root");
  }

  //ofstream frunlumi;
  //frunlumi.open("run_lumi_events_runC.txt");

  double ntot = 0;
  double nPhoPrompt = 0;
  double ntotLep = 0;
  double nolep = 0;
  double totpos = 0;
  double totneg = 0;

  //showershape correction
  //TFile* fss_rw = TFile::Open("external/transformation_Moriond17_AfterPreApr_v1.root");
  TFile* fss_rw = TFile::Open("external/transformation_pho_presel_BDTUpto6000.root");
  TGraph *tgr[14];
  tgr[0] = (TGraph*) fss_rw->Get("transfEtaWidthEB");
  tgr[1] = (TGraph*) fss_rw->Get("transfS4EB");
  tgr[2] = (TGraph*) fss_rw->Get("transffull5x5R9EB");
  tgr[3] = (TGraph*) fss_rw->Get("transffull5x5sieieEB");
  tgr[4] = (TGraph*) fss_rw->Get("transffull5x5sieipEB");
  tgr[5] = (TGraph*) fss_rw->Get("transfrhoEB");
  tgr[6] = (TGraph*) fss_rw->Get("transfPhiWidthEB");

  tgr[7] = (TGraph*) fss_rw->Get("transfEtaWidthEE");
  tgr[8] = (TGraph*) fss_rw->Get("transfS4EE");
  tgr[9] = (TGraph*) fss_rw->Get("transffull5x5R9EE");
  tgr[10] = (TGraph*) fss_rw->Get("transffull5x5sieieEE");
  tgr[11] = (TGraph*) fss_rw->Get("transffull5x5sieipEE");
  tgr[12] = (TGraph*) fss_rw->Get("transfrhoEE");
  tgr[13] = (TGraph*) fss_rw->Get("transfPhiWidthEE");

  for (Long64_t ev = ev1; ev < ev2; ++ev) {
    //  for (Long64_t ev = 0; ev < 50000; ++ev) {
    if (ev % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data1.GetEntriesFast());
    
    //cout << ev << endl;

    data1.GetEntry(ev);
    readggtree(data1);

    if (!isPVGood) continue; //select good vertex

    //vertex
    vz_ = vz;
    vx_ = vx;
    vy_ = vy;
    nVtx_ = nVtx;
    run_ = run;
    event_ = event;
    lumi_ = lumis;
    rho_ = rho;

    /*
    pdfWeight_ = pdfWeight;

    //weight for NNPDF_nlo_nf5_pdfas
    for (int i = 0; i < 102; i++) {
      nnpdf_nlo_wei_[i] = pdfSystWeight[i+1];
    }

    nnpdf_nlo_scale_wei_[0] = pdfSystWeight[104];
    nnpdf_nlo_scale_wei_[1] = pdfSystWeight[105];
    nnpdf_nlo_scale_wei_[2] = pdfSystWeight[106];
    nnpdf_nlo_scale_wei_[3] = pdfSystWeight[107];
    nnpdf_nlo_scale_wei_[4] = pdfSystWeight[108];
    nnpdf_nlo_scale_wei_[5] = pdfSystWeight[110];
    nnpdf_nlo_scale_wei_[6] = pdfSystWeight[112];

    //weight for PDF4LHC15_nlo_100_pdfas
    for (unsigned int i = 0; i <= 102; i++) {
      PDF4LHC_wei_[i] = pdfSystWeight[i+465];
    }
    */
    //PDF4LHC_scale_wei_[0] = pdfSystWeight[105];
    //PDF4LHC_scale_wei_[1] = pdfSystWeight[106];
    //PDF4LHC_scale_wei_[2] = pdfSystWeight[107];
    //PDF4LHC_scale_wei_[3] = pdfSystWeight[108];
    //PDF4LHC_scale_wei_[4] = pdfSystWeight[110];
    //PDF4LHC_scale_wei_[5] = pdfSystWeight[112];
    

    float lheElePhi1 = -99., lheElePhi2 = -99.;
    float lheMuPhi1 = -99., lheMuPhi2 = -99.;

    TLorentzVector  lheEle_Vec1, lheEle_Vec2;
    TLorentzVector  lheMu_Vec1, lheMu_Vec2;
    TLorentzVector  lhePho_Vec, lheZ_Vec, lheX_Vec;

    for (int i = 0; i < 3; i++) {
      lheJetPt[i] = -99.;
      lheJetEta[i] = -99.;
    }

    nlheJet = 0;
    //get lhe pt for lepton
    if (data1.HasMC()) {
      lhePho1 = genPho1;
      lhePhoEta1 = genPhoEta1;

      dRlhePhoEle1 = deltaR(genEleEta1, genElePhi1, genPhoEta1, genPhoPhi1);
      dRlhePhoEle2 = deltaR(genEleEta2, genElePhi2, genPhoEta1, genPhoPhi1);
      dRlhePhoMu1 = deltaR(genMuEta1, genMuPhi1, genPhoEta1, genPhoPhi1);
      dRlhePhoMu2 = deltaR(genMuEta2, genMuPhi2, genPhoEta1, genPhoPhi1);


      lheEle_Vec1.SetPtEtaPhiM(genEle1, genEleEta1, genElePhi1, 0.511*0.001);
      lheEle_Vec2.SetPtEtaPhiM(genEle2, genEleEta2, genElePhi2, 0.511*0.001);
      
      lheMu_Vec1.SetPtEtaPhiM(genMu1, genMuEta1, genMuPhi1, 0.105);
      lheMu_Vec2.SetPtEtaPhiM(genMu2, genMuEta2, genMuPhi2, 0.105);

      lhePho_Vec.SetPtEtaPhiM(genPho1, genPhoEta1, genPhoPhi1, 0);

      for (int i =0 ; i < nLHEJet; i++) {
	if (lheJetPt_[i] < 30.) continue;
	if ( abs(lheJetEta_[i]) > 2.4) continue;
	if (deltaR(genPhoEta1, genPhoPhi1,lheJetEta_[i],lheJetPhi_[i]) < 0.4) continue;
	if (genEle1 > 0. && genEle2 >0.) {
	  if (deltaR(genEleEta1, genElePhi1,lheJetEta_[i],lheJetPhi_[i]) < 0.4) continue;
	  if (deltaR(genEleEta2, genElePhi2,lheJetEta_[i],lheJetPhi_[i]) < 0.4) continue;
	}
	if (genMu1 > 0. && genMu2 >0.) {
	  if (deltaR(genEleEta1, genElePhi1,lheJetEta_[i],lheJetPhi_[i]) < 0.4) continue;
	  if (deltaR(genEleEta2, genElePhi2,lheJetEta_[i],lheJetPhi_[i]) < 0.4) continue;
	}
	  

	lheJetPt[nlheJet] = lheJetPt_[i];
	lheJetEta[nlheJet] = lheJetEta_[i];
	//if (nlheJet > 0) cout << "lheJet Pt: " << lheJetPt_[i] << "\t" << lheJetPt[nlheJet] << endl;
	nlheJet++;

      }

      
      if (genEle1 > genEle2) {
	lheEle1 = genEle1;
	lheEleEta1 = genEleEta1;
	lheElePhi1 = genElePhi1;

	lheEle2 = genEle2;
	lheEleEta2 = genEleEta2;
	lheElePhi2 = genElePhi2;

      }
      else {
	lheEle1 = genEle2;
	lheEleEta1 = genEleEta2;
        lheElePhi1 = genElePhi2;

	lheEle2 = genEle1;
	lheEleEta2 = genEleEta1;
	lheElePhi2 = genElePhi1;

      }
      
      if (genMu1 > genMu2) {
	lheMu1 = genMu1;
	lheMuEta1 = genMuEta1;
        lheMuPhi1 = genMuPhi1;

	lheMu2 = genMu2;
	lheMuEta2 = genMuEta2;
	lheMuPhi2 = genMuPhi2;

      }
      else {
	lheMu1 = genMu2;
	lheMuEta1 = genMuEta2;
        lheMuPhi1 = genMuPhi2;

	lheMu2 = genMu1;
	lheMuEta2 = genMuEta1;
        lheMuPhi2 = genMuPhi1;
      }
    }

    //get lhe Z and Zg
    if ( lheEle1 > 0. && lheEle2 > 0.) {
      lheZ_Vec = lheEle_Vec1 + lheEle_Vec2;
      lhe_Zm = lheZ_Vec.M();
      lhe_type = 11;
    }
    else if (  lheMu1 > 0. && lheMu2 > 0.) {
      lheZ_Vec = lheMu_Vec1 + lheMu_Vec2;
      lhe_Zm = lheZ_Vec.M();
      lhe_type = 13;
    }
    else {
      lhe_Zm = -99.;
      lhe_type = -1;
    }

    if (lhePho1 > 0. && lhe_Zm > 0.) {

      lheX_Vec = lheZ_Vec + lhePho_Vec;
      lhe_Xm = lheX_Vec.M();
    }
    else {
      lhe_Xm = -99.;
    }

    
    //PU reweighting for MC
    if (data1.HasMC()) {
      float* puTrue = data1.GetPtrFloat("puTrue");
      puweigj = (float) puCalcGJ_69nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_69p2nb = (float) puCalcGJ_69p2nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_65nb = (float) puCalcGJ_65nb.GetWeight(run, puTrue[1]); // in-time PU
      puweigj_63nb = (float) puCalcGJ_63nb.GetWeight(run, puTrue[1]); // in-time PU
      //puweisj = (float) puCalcSJ.GetWeight(run, puTrue[1]); // in-time PU
      
      genWeight_ = (genWeight > 0) ? 1. : -1.;
      ntot += genWeight_;
    }

    //generated level
    bool skipEvent = false;


    if (data1.HasMC()) {
    int genpho = 0;
    TLorentzVector genEle[5], genMu[5], GenPho;
    TLorentzVector genlep[5];
    TLorentzVector genfsr_lep1[10], genfsr_lep2[10];
    TLorentzVector genfsr_ele1[10], genfsr_ele2[10];
    TLorentzVector genfsr_mu1[10], genfsr_mu2[10];
    TLorentzVector genZ, genX, genZll;
    
    ngenEle = 0;
    ngenMu = 0;

    lep1fsr = 0;
    lep2fsr = 0;
    ngenlep = 0;
    genPhoEt = -99.;
    genPhoEta = -99.;
    genPhoY = -99.;
    genPhoPhi = -99.;
    gendRPhoLep1 = -99.;
    gendRPhoLep2 = -99.;
    genCalIso03 = -99.;
    genCalIso04 = -99.;
    genZm = -99.;
    mcZm = -99.;
    no_fsr_pho = 0;
    for (int i = 0; i < 50; i++) {
      fsr_gen_pho_pt[i] = -99.;
      fsr_pho_lep_dr1[i] = -99.;
      fsr_pho_lep_dr2[i] = -99.;
    }
	
    for (int i = 0; i < nMC; i++) {
      if ( abs(mcPID[i])==11 && ((mcStatusFlag[i]>>0)&1)==1 ) {
	//if ( abs(mcPID[i])==11 && ((mcStatusFlag[i]>>8)&1)==1 ) { //before FSR
	genlepPt[ngenlep] = mcPt[i];
	genlepEta[ngenlep] = mcEta[i];
	genlepPhi[ngenlep] = mcPhi[i];
	genlepCh[ngenlep] = (mcPID[i]==11) ? -1: 1;
	genlep[ngenEle].SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.511*0.001);
	//genZll.SetPtEtaPhiM(mcMomPt[i], mcMomEta[i], mcMomPhi[i], mcMomMass[i]);
	//mcZm = mcMomMass[i];
	ngenEle++;
	ngenlep++;
      }

      if (abs(mcPID[i])==13 && ((mcStatusFlag[i]>>0)&1)==1) {
	//if (abs(mcPID[i])==13 && ((mcStatusFlag[i]>>8)&1)==1) { //before FSR
	genlepPt[ngenlep] = mcPt[i];
	genlepEta[ngenlep] = mcEta[i];
	genlepPhi[ngenlep] = mcPhi[i];
	genlepCh[ngenlep] = (mcPID[i]==13) ? -1: 1;
	genlep[ngenMu].SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.105);
	//genZll.SetPtEtaPhiM(mcMomPt[i], mcMomEta[i], mcMomPhi[i], mcMomMass[i]);
	//mcZm = mcMomMass[i];
	ngenMu++;
	ngenlep++;
	}
    
      //FSR
      //if ( mcPID[i]==22 && (abs(mcMomPID[i])==11 ||abs(mcMomPID[i])==13) ) {
      if ( mcPID[i]==22) {
	//cout << "calculate dr1 and dr2 " << endl;
	float dr1 = deltaR(genlep[0].Eta(), genlep[0].Phi(), mcEta[i], mcPhi[i]);
	float dr2 = deltaR(genlep[1].Eta(), genlep[1].Phi(), mcEta[i], mcPhi[i]);
	if (dr1 < 0.1) { 
	  genfsr_lep1[lep1fsr].SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.);
	  lep1fsr++;
	}
	if ( dr2 < 0.1) {
	  genfsr_lep2[lep1fsr].SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.);
	  lep2fsr++;
	}
	}

    } //end loop on nMC


      //put pack FSR photon to lep
    if (ngenEle >= 1 || ngenMu >= 1) {
      ntotLep++;
      if ( lep1fsr >= 1) {
	for ( int ifsr = 0; ifsr < lep1fsr; ifsr++) 
	  genlep[0] = genlep[0] + genfsr_lep1[ifsr];
      }
      if ( lep2fsr >= 1) {
	for ( int ifsr = 0;ifsr < lep2fsr; ifsr++)
	  genlep[1] = genlep[1] + genfsr_lep2[ifsr];
      }
    }
    

    //order in lepton pt descending
    float tmp_pt, tmp_eta;
    if (ngenEle > 1 || ngenMu > 1) {
      if (genlep[0].Pt() < genlep[1].Pt()) {
	tmp_pt = genlep[1].Pt();
	genlepPt[1] = genlep[0].Pt();
	genlepPt[0] = tmp_pt;

	tmp_eta = genlep[1].Eta();
        genlepEta[1] = genlep[0].Eta();
        genlepEta[0] = tmp_eta;

      }
      else {
	genlepPt[0] = genlep[0].Pt();
	genlepPt[1] = genlep[1].Pt();
	genlepEta[0]= genlep[0].Eta();
	genlepEta[1]= genlep[1].Eta();
      }
    }

      
    genZtype = 0;
    //cout << ngenEle << "\t " << ngenMu << endl;
    if (ngenEle > 1 || ngenMu > 1) {
      genZ = genlep[0]+genlep[1];
      genZm = genZ.M();
      genZpt = genZ.Pt();
      genZy = genZ.Rapidity();
      genZch = genlepCh[0]+genlepCh[1];
    }
    if (ngenEle > 1)
      genZtype = 11; 
    else if ( ngenMu > 1)
      genZtype = 13;
    else {
      genZm = -999.;
      genZpt = -999.;
      genZy = -999.;
      genZch = -99;
    }
    

    //selection gen photon
    genpho = 0;
    vector<float> genphoet;
    vector<int> genindex;

    for (int i = 0; i < nMC; i++) {
      //if (mcPID[i]==22 && ( ((mcStatusFlag[i]>>0)&1)==1 || ((mcStatusFlag[i]>>1)&1) == 1) && (abs(mcMomPID[i])==11 || abs(mcMomPID[i])==13 || abs(mcMomPID[i])==15)) {
      if (mcPID[i]==22 && ( ((mcStatusFlag[i]>>0)&1)==1 || ((mcStatusFlag[i]>>1)&1) == 1) ) {
	//if (mcPID[i]==22 && ((mcStatusFlag[i]>>0)&1)==1) {
	nPhoPrompt++;
	//if (mcCalIsoDR04[i] > 5.) skipEvent = true;
	if (mcCalIsoDR04[i] > 5.) continue;
	if ( mcPt[i] < 15.) continue;
	//if ( abs(mcEta[i]) > 2.5) continue;
	float dr1 = deltaR(genlep[0].Eta(), genlep[0].Phi(), mcEta[i], mcPhi[i]);
	float dr2 = deltaR(genlep[1].Eta(), genlep[1].Phi(), mcEta[i], mcPhi[i]);
	if ( dr1 < 0.7 || dr2 < 0.7 ) continue;

	genphoet.push_back(mcPt[i]);
	genindex.push_back(i);
      }
    }


    //choose highest photon pt
    float maxtmp = 0.;
    int igenpho = -1;
    for (vector<int>::iterator igen = genindex.begin(); igen != genindex.end(); igen++) {
      if (maxtmp < mcPt[*igen]) {
	  maxtmp = mcPt[*igen];
	  igenpho = *igen;
	  genpho++;
      }
    }

    if (igenpho >= 0) {
      genPhoEt = mcPt[igenpho];
      genPhoEta = mcEta[igenpho];
      genPhoPhi = mcPhi[igenpho];
      genCalIso03 = mcCalIsoDR03[igenpho];
      genCalIso04 = mcCalIsoDR04[igenpho];
      GenPho.SetPtEtaPhiM(mcPt[igenpho], mcEta[igenpho], mcPhi[igenpho], 0.);
      genPhoY = GenPho.Rapidity();
      gendRPhoLep1 = deltaR(genlep[0].Eta(), genlep[0].Phi(), GenPho.Eta(), GenPho.Phi());
      gendRPhoLep2 = deltaR(genlep[1].Eta(), genlep[1].Phi(), GenPho.Eta(), GenPho.Phi());
    }


    if (genZm>0 && genPhoEt>0.) {
      genX = genZ + GenPho;
      //genX = genZll + GenPho;
      genXm = genX.M();
      genXpt = genX.Pt();
      if (genWeight > 0) totpos++;
      else totneg++;
    }
    else {
      genXm = -99.;
      genXpt = -99.;
    }

    if (inputfile.Contains("_ZZTo4L.") || inputfile.Contains("_WZTo3LNu.") ) {
      if (ngenEle >= 3 || ngenMu >= 3) {
	genX = genlep[0] + genlep[1] + genlep[2];
	genXm = genX.M();
	genXpt = genX.Pt();
      }
      else {
	genXm = -99.;
	genXpt = -99.;
      }
    }


    //select gen jet
    int nmaxjet = 60;
    for (int i = 0; i < nmaxjet; i++) {
      genJetPt[i] = -99.;
      genJetEta[i] = -99.;
    }
    
    ngenjet_ = 0;
    for (int i = 0; i < ngenjet; i++) {
      if ( genjetPt_all[i] < 30.) continue;
      if ( abs(genjetEta_all[i]) > 2.4) continue;
      if (genlep[0].Pt() > 0 && deltaR(genlep[0].Eta(), genlep[0].Phi(), genjetEta_all[i], genjetPhi_all[i]) < 0.4) continue;
      if (genlep[1].Pt() > 0 && deltaR(genlep[1].Eta(), genlep[1].Phi(), genjetEta_all[i], genjetPhi_all[i]) < 0.4) continue;
      if (genPhoEt>0. && deltaR(genPhoEta, genPhoPhi, genjetEta_all[i], genjetPhi_all[i]) < 0.4) continue;

      genJetPt[ngenjet_] = genjetPt_all[i];
      genJetEta[ngenjet_] = genjetEta_all[i];
      ngenjet_++;
    }
    
    outtreeGen->Fill();
    }
    
    //-------- done gen level ---------//

    //------------------------------------------//
    //--------------- reco levl ---------------- //
    
    //initial value for output variables
    leptType        = -99;
    weight          = -99;
    aMCNLO          = 0;
    nPU             = -99;
    //nVert           = -99;
    lept0_pt        = -99.;
    lept0_eta       = -99.;
    lept0_phi       = -99.;
    lept1_sceta     = -99.;
    lept0_miniRelIso = -99.;
    lept0_trkIso    = -99.;
    lept0_pdgId     = -99;
    lept0_normalizedChi2 = -99.;
    lept0_nValidMuonHits = -99.;
    lept0_nMatchedStations = -99;
    lept0_nValidPixelHits = -99;
    lept0_nTrackerLayers = -99;
    lept0_muonBestTrack_dxyVTX = -99.;
    lept0_muonBestTrack_dzVTX = -99.;
    lept0_ptError   = -99.;
    lept0_sigmaIetaIeta = -99.;
    lept0_dEtaIn    = -99.;
    lept0_dPhiIn    = -99.;
    lept0_hOverE    = -99.;
    lept0_ooEmooP   = -99.;
    lept0_d0        = -99.;
    lept0_dz        = -99.;
    lept0_mva = -99.;
    lept0_chiso = -99.;
    lept0_phoiso = -99.;
    lept0_neuiso = -99.;
    lept0_sigEOverE = -99.;
    lept0_EnErr     = -99.;
    lept0_En        = -99.;
    lept0_expectedMissingInnerHits = -99;
    lept0_RecoSF = 1.;
    lept0_SelSF = 1.;
    lept0_trigSF = 1.;

    lept0_RecoSF_Err = 0;
    lept0_SelSF_Err = 0;
    lept0_trigSF_Err = 0;

    lept1_pt        = -99.;
    lept1_eta       = -99.;
    lept1_phi       = -99.;
    lept1_sceta     = -99.;
    lept1_miniRelIso = -99.;
    lept1_trkIso    = -99.;
    lept1_pdgId     = -99;
    lept1_normalizedChi2 = -99.;
    lept1_nValidMuonHits = -99;
    lept1_nMatchedStations = -99;
    lept1_nValidPixelHits = -99;
    lept1_nTrackerLayers = -99;
    lept1_muonBestTrack_dxyVTX = -99.;
    lept1_muonBestTrack_dzVTX = -99.;
    lept1_ptError   = -99.;
    lept1_sigmaIetaIeta = -99.;
    lept1_dEtaIn    = -99.;
    lept1_dPhiIn    = -99.;
    lept1_hOverE    = -99.;
    lept1_ooEmooP   = -99.;
    lept1_d0        = -99.;
    lept1_dz        = -99.;
    lept1_expectedMissingInnerHits = -99;
    lept1_EnErr     = -99.;
    lept1_mva = -99.;
    lept1_chiso = -99.;
    lept1_phoiso = -99.;
    lept1_neuiso = -99.;
    lept1_sigEOverE = -99.;
    lept1_En        = -99.;
    lept1_RecoSF = 1.;
    lept1_SelSF = 1.;
    lept1_trigSF = 1.;
    lept_dzSF = 1.;
    lept_dzSF_Err = 0.;
    lept1_RecoSF_Err = 0.;
    lept1_SelSF_Err = 0.;
    lept1_trigSF_Err = 0.;
    deltaR_lept     = -99.;
    deltaPhi_lept   = -99.;
    isEB            = 0;
    isEE            = 0;
    gamma_SF        = 1.;
    gamma_SF_Err    = 0.;
    gamma_CSEV_SF        = 1.;
    gamma_CSEV_SF_Err    = 0.;
    gamma_pt        = -99.;
    gamma_eta       = -99.;
    gamma_sceta  = -99.;
    gamma_phi       = -99.;
    gamma_HoverE    = -99.;
    gamma_sigmaIetaIeta = -99.;
    gamma_mva       = -99.;
    gamma_ssmva     = -99.;
    gamma_ssmva_norw = -99.;
    gamma_R9        = -99.;
    gamma_ChIso     = -99.;
    gamma_NeuIso     = -99.;
    gamma_PhoIso     = -99.;
    gamma_ChWorstIso= -99.;
    gamma_PixSeed   = -99.;
    gamma_EleVeto   = -99.;
    gamma_En        = -99.;
    gamma_EnErr     = -99.;
    gamma_geniso    = -99.;
    gamma_SCEtaWidth = -99.;
    gamma_s4Full5x5 = -99.;
    gamma_SCEtaWidth_rw = -99.;
    gamma_s4Full5x5_rw = -99.;
    gamma_R9_rw        = -99.;
    gamma_sigmaIetaIeta_rw = -99.;
    deltaR_PhoLept0 = -99.;
    deltaR_PhoLept1 = -99.;
    deltaPhi_PhoLept0 = -99.;
    deltaPhi_PhoLept1 = -99.;
    z_pt            = -99.;
    z_eta           = -99.;
    z_phi           = -99.;
    z_y = -99.;
    z_mass          = -99.;
    boss_pt         = -99.;
    boss_eta        = -99.;
    boss_phi        = -99.;
    boss_mass       = -99.;
    boss_massrel    = -99.;


    //trigger
    trig_Ele27_WPTight = hlt>>4&1;
    trig_Ele23_Ele12 = ( (hlt>>40&1) ==1 || (hlt>>5&1) == 1) ? 1 : 0;
    //if ( ((hlt>>40)&1) || ((hlt>>5)&1) ) trig_Ele23_Ele12 = 1;
    //else trig_Ele23_Ele12 = 0;
    trig_Ele23_Ele12_nonDZ = (hlt>>40)&1;
    trig_Ele23_Ele12_DZ = (hlt>>5)&1;

    
    trig_DoubleEle33 = hlt>>11&1;

    trig_Mu50 = hlt >> 21&1;
    if ( (hlt>>14&1) || (hlt>>15&1) || (hlt>>41&1) || (hlt>>42&1))
      trig_Mu17_Mu8 = 1;
    else trig_Mu17_Mu8 = 0;

    trig_Mu17_Mu8_nonDZ = (hlt>>41&1) || (hlt>>42&1);
    trig_Mu17_Mu8_DZ = (hlt>>14&1) || (hlt>>15&1);

    trig_doublePho60 = hltPho >> 22&1;
    //if ( !trig_Ele23_Ele12 && !trig_Mu17_Mu8) continue;

    TLorentzVector ele[5], mu[5], pair;
    int ne = 0;
    int nmu = 0;
    int npho = 0;

    float eleRecoSF_[5];
    float eleSF_[5];
    float err_eleRecoSF_[5];
    float err_eleSF_[5];

    vector<int> eleIndex;
    vector<int> selectedEle;
    //Electron_SMZg(selectedEle);
    passEleId(selectedEle);
    //NoEleId(selectedEle);
    //cout << "passing electron tight id" << endl;

    if ( selectedEle.size() > 1 ) {

      int match_wwqq = 0;
      //match 1lep to genlep, 1lep to jet
      if ( inputfile.Contains("WWToLNuQQ") ) {
	int matchlepjet[5];
	for (vector<int>::const_iterator ie = selectedEle.begin(); ie != selectedEle.end(); ie++) {
	  matchlepjet[ne] = eleMatcher(data1, *ie);
	  ne++;
	}

	if ( (matchlepjet[0]==0 &&  matchlepjet[1]==1) || (matchlepjet[0]== 1 &&  matchlepjet[1]==0) )
	  //cout << "matching for WWToLNuQQ sample" << endl;
	  match_wwqq = 1;
      } //end of matching for WWToLNuQQ
      ne = 0;

    for (vector<int>::const_iterator ie = selectedEle.begin(); ie != selectedEle.end(); ie++) {
      //cout << "check mc or not" << endl;
      if (mc) {
	if (inputfile.Contains("WWToLNuQQ")) {
	  //cout << "matching for WWToLNuQQ sample: " << match_wwqq << endl;
	    if (match_wwqq != 1) continue;
	  }
	else {
	  bool match = eleMatcher(data1, *ie);
	  if ( !match) continue;
	}
      }

      //TRandom3 *gen = new TRandom3(0);
      //elePt[*ie] *= (gen->Gaus(1,eleResol_rho_up[*ie]));
      //elePt[*ie] *= eleScale_gain_dn[*ie];

      eleIndex.push_back(*ie);
      ele[ne].SetPtEtaPhiM(elePt[*ie],eleEta[*ie],elePhi[*ie], 0.511*0.001);

      ne++;
    }
    }  //end of electron loop

    //muon loop
    int indexMu1 = -99, indexMu2 = -99;
    vector<int> muIndex;
    nmu = 0;

    vector<int> selectedMu;
    passMuonId(selectedMu);

    //cout << "-----muon loop---" << endl;

    if ( selectedMu.size() > 1 ) {

      int match_wwqq = 0;
      //match 1lep to genlep, 1lep to jet
      if ( inputfile.Contains("WWToLNuQQ") ) {
	int matchlepjet[5];
	for (vector<int>::const_iterator imu = selectedMu.begin(); imu != selectedMu.end(); imu++) {
	  matchlepjet[nmu] = muonMatcher(data1, *imu);
	  nmu++;
	}

	if ( (matchlepjet[0]==0 &&  matchlepjet[1]==1) || (matchlepjet[0]== 1 &&  matchlepjet[1]==0) )
	  match_wwqq = 1;
      } //end of matching for WWToLNuQQ
      nmu = 0;

    for (vector<int>::const_iterator imu = selectedMu.begin(); imu != selectedMu.end(); imu++) {
      //cout << "check mc or not" << endl;
      if (mc) {
	if (inputfile.Contains("WWToLNuQQ")) {
	  if (match_wwqq !=1) continue;
	}
	else {
	  bool match = muonMatcher(data1, *imu);
	  if ( !match) continue;
	}
      }

      //muon rochester correction
      double SF      = 0.;
      float  muCalibPt      = 0.;

      //for mc
      if (mc) {
	float genPt = -99.;
        for (Int_t j=0; j<nMC; ++j) {
          if (abs(mcPID[j]) != 13) continue;
          if (deltaR(mcEta[j], mcPhi[j], muEta[*imu], muPhi[*imu]) < 0.1) genPt = mcPt[j];
        }

	int nl = muTrkLayers[*imu];
	float rnd1 = gRandom->Rndm();
	if (genPt > 0) {
          SF      = rc.kScaleFromGenMC(muCh[*imu], muPt[*imu], muEta[*imu], muPhi[*imu], genPt, 0, 0);
	}
	else {
	  SF      = rc.kScaleAndSmearMC(muCh[*imu], muPt[*imu], muEta[*imu], muPhi[*imu], nl, rnd1, 0, 0);
	}
      }
      else {
	SF = rc.kScaleDT(muCh[*imu], muPt[*imu], muEta[*imu], muPhi[*imu], 0, 0); //for data
      }
	
      muCalibPt =  muPt[*imu]*SF;
      
      muIndex.push_back(*imu);
      mu[nmu].SetPtEtaPhiM(muCalibPt, muEta[*imu], muPhi[*imu], 0.105);
      nmu++;
    }
    }
    
    //cout << "done muon loop" << endl;

    //cout << "nmu = " << nmu << "\t ne: " << ne << endl;

    if (nmu <=1 && ne <=1) continue;  //event has at least 2e or 2m

    if (ne > 1) {
      if ( elePt[eleIndex[0]] < 2.) continue; //for SM-Zg
      if ( elePt[eleIndex[1]] < 2.) continue; //for SM-Zg

      //cout << "ele variable: " << endl;
      leptType = 11;
      lept0_pt = elePt[eleIndex[0]];
      lept0_eta = eleEta[eleIndex[0]];
      lept0_sceta = eleSCEta[eleIndex[0]];
      lept0_phi = elePhi[eleIndex[0]];
      //lept0_miniRelIso = elePFMiniIso[eleIndex[0]];
      lept0_pdgId = (eleCharge[eleIndex[0]]==1) ? -11 : 11;
      //lept0_normalizedChi2 = eleGSFChi2[eleIndex[0]];
      lept0_sigmaIetaIeta = eleSigmaIEtaIEta_Full5x5[eleIndex[0]];
      lept0_dEtaIn = eledEtaAtVtx[eleIndex[0]];
      lept0_dPhiIn = eledPhiAtVtx[eleIndex[0]];
      lept0_hOverE = eleHoverE[eleIndex[0]];
      lept0_ooEmooP = eleEoverPInv[eleIndex[0]];
      lept0_d0 = eleD0[eleIndex[0]];
      lept0_dz = eleDz[eleIndex[0]];
      //lept0_mva = eleIDMVAHZZ[eleIndex[0]];
      lept0_chiso = elePFChIso[eleIndex[0]];
      lept0_phoiso = elePFPhoIso[eleIndex[0]];
      lept0_neuiso = elePFNeuIso[eleIndex[0]];
      //lept0_SIP = eleSIP[eleIndex[0]];
      //lept0_expectedMissingInnerHits = eleMissHits[eleIndex[0]];
      lept0_En = eleEn[0];

      //cout << "scale infor" << endl;
      eleScale_stat_up_1 = eleScale_stat_up[eleIndex[0]];
      eleScale_stat_dn_1 = eleScale_stat_dn[eleIndex[0]];
      eleScale_syst_up_1 = eleScale_syst_up[eleIndex[0]];
      eleScale_syst_dn_1 = eleScale_syst_dn[eleIndex[0]];
      eleScale_gain_up_1 = eleScale_gain_up[eleIndex[0]];
      eleScale_gain_dn_1 = eleScale_gain_dn[eleIndex[0]];
      eleResol_rho_up_1 = eleResol_rho_up[eleIndex[0]];
      eleResol_rho_dn_1 = eleResol_rho_dn[eleIndex[0]];
      eleResol_phi_up_1 = eleResol_phi_up[eleIndex[0]];
      eleResol_phi_dn_1 = eleResol_phi_dn[eleIndex[0]];


      //cout  << "trailing ele" << endl;
      //variable for trailing lepton
      lept1_pt = elePt[eleIndex[1]];
      lept1_eta = eleEta[eleIndex[1]];
      lept1_sceta = eleSCEta[eleIndex[1]];
      lept1_phi = elePhi[eleIndex[1]];
      //lept1_miniRelIso = elePFMiniIso[eleIndex[1]];
      lept1_pdgId = (eleCharge[eleIndex[1]]==1) ? -11 : 11;
      //lept1_normalizedChi2 = eleGSFChi2[eleIndex[1]];
      lept1_sigmaIetaIeta = eleSigmaIEtaIEta_Full5x5[eleIndex[1]];
      lept1_dEtaIn = eledEtaAtVtx[eleIndex[1]];
      lept1_dPhiIn = eledPhiAtVtx[eleIndex[1]];
      lept1_hOverE = eleHoverE[eleIndex[1]];
      lept1_ooEmooP = eleEoverPInv[eleIndex[1]];
      lept1_d0 = eleD0[eleIndex[1]];
      lept1_dz = eleDz[eleIndex[1]];
      lept1_expectedMissingInnerHits = eleMissHits[eleIndex[1]];
      //lept1_mva = eleIDMVAHZZ[eleIndex[1]];
      lept1_chiso = elePFChIso[eleIndex[1]];
      lept1_phoiso = elePFPhoIso[eleIndex[1]];
      lept1_neuiso = elePFNeuIso[eleIndex[1]];
      //lept1_SIP = eleSIP[eleIndex[1]];
      //lept1_EnErr = eleEcalEnErr[1];
      lept1_En = eleEn[1];

      //cout << "get sf for ele " << endl;
      if (mc) {
	lept0_RecoSF = ele_RecoSF(lept0_pt, lept0_sceta);
	lept0_SelSF = ele_SelSF(lept0_pt, lept0_sceta);
	lept0_trigSF = ele_TrigSF_leg1(lept0_pt, lept0_sceta);;
	
	lept0_RecoSF_Err = ele_RecoSF_err(lept0_pt, lept0_sceta);
	lept0_SelSF_Err = ele_SelSF_err(lept0_pt, lept0_sceta);
	lept0_trigSF_Err = ele_TrigSF_leg1_err(lept0_pt, lept0_sceta);

	lept1_RecoSF = ele_RecoSF(lept1_pt, lept1_sceta);
	lept1_SelSF = ele_SelSF(lept1_pt, lept1_sceta);
	lept1_trigSF = ele_TrigSF_leg2(lept1_pt, lept1_sceta);;
	
	lept1_RecoSF_Err = ele_RecoSF_err(lept1_pt, lept1_sceta);
	lept1_SelSF_Err = ele_SelSF_err(lept1_pt, lept1_sceta);
	lept1_trigSF_Err = ele_TrigSF_leg2_err(lept1_pt, lept1_sceta);

	lept_dzSF = ele_trigSF_DZ(abs(lept0_eta), abs(lept1_eta));
	lept_dzSF_Err = ele_trigSF_DZ_err(abs(lept0_eta), abs(lept1_eta));

	//cout << "Err of Sel SF for ele1: " << lept0_SelSF_Err << "\t ele2: " << lept1_SelSF_Err << endl;
	  
      }

      eleScale_stat_up_2 = eleScale_stat_up[eleIndex[1]];
      eleScale_stat_dn_2 = eleScale_stat_dn[eleIndex[1]];
      eleScale_syst_up_2 = eleScale_syst_up[eleIndex[1]];
      eleScale_syst_dn_2 = eleScale_syst_dn[eleIndex[1]];
      eleScale_gain_up_2 = eleScale_gain_up[eleIndex[1]];
      eleScale_gain_dn_2 = eleScale_gain_dn[eleIndex[1]];
      eleResol_rho_up_2 = eleResol_rho_up[eleIndex[1]];
      eleResol_rho_dn_2 = eleResol_rho_dn[eleIndex[1]];
      eleResol_phi_up_2 = eleResol_phi_up[eleIndex[1]];
      eleResol_phi_dn_2 = eleResol_phi_dn[eleIndex[1]];


      pair = ele[0] + ele[1];
      z_mass = pair.M();
      z_pt = pair.Pt();
      z_y = pair.Rapidity();
      z_eta = pair.Eta();
      z_phi = pair.Phi();
      z_charge = eleCharge[eleIndex[0]] + eleCharge[eleIndex[1]];

      deltaPhi_lept = deltaPhi(ele[0].Phi(),ele[1].Phi());
      deltaR_lept = deltaR(ele[0].Eta(), ele[0].Phi(),ele[1].Eta(), ele[1].Phi());

      if (pair.M() > 60 && pair.M() < 120) {
	if (z_charge == 0) heemass->Fill(pair.M());
	else heemassSS->Fill(pair.M());
      }

    }

    else if (nmu > 1) {
      indexMu1 = muIndex[0];
      indexMu2 = muIndex[1];
      if (muPt[indexMu1] < 2.) continue;
      if (muPt[indexMu2] < 2.) continue;

      //cout << "get infor of muon variable" << endl;
      leptType = 13;
      lept0_pt = muPt[indexMu1];
      lept0_eta = muEta[indexMu1];
      lept0_phi = muPhi[indexMu1];
      lept0_pdgId = (muCh[indexMu1]==1) ? -13 : 13;
      /*
      lept0_miniRelIso = muPFMiniIso[indexMu1];
      lept0_d0 = muD0[indexMu1];
      lept0_dz = muDz[indexMu1];
      lept0_trkIso = muIsoTrk[indexMu1];
      lept0_pdgId = (muCh[indexMu1]==1) ? -13 : 13;
      lept0_normalizedChi2 = muChi2NDF[indexMu1];
      lept0_nValidMuonHits = muMuonHits[indexMu1];
      lept0_nMatchedStations = muStations[indexMu1];
      lept0_nValidPixelHits = muPixelHits[indexMu1];
      lept0_nTrackerLayers = muTrkLayers[indexMu1];
      lept0_muonBestTrack_dxyVTX = muD0[indexMu1];
      lept0_muonBestTrack_dzVTX = muDz[indexMu1];
      lept0_ptError = muBestTrkPtError[indexMu1];
      lept0_chiso = muPFChIso03[indexMu1];
      lept0_phoiso = muPFPhoIso03[indexMu1];
      lept0_neuiso = muPFNeuIso03[indexMu1];
      lept0_SIP = muSIP[indexMu1];
      lept0_sigEOverE= muBestTrkPtError[indexMu1]/muBestTrkPt[indexMu1];
      */
      
      lept1_pt = muPt[indexMu2];
      lept1_eta = muEta[indexMu2];
      lept1_phi = muPhi[indexMu2];
      /*
      lept1_miniRelIso = muPFMiniIso[indexMu2];
      lept1_trkIso = muIsoTrk[indexMu2];
      lept1_pdgId = (muCh[indexMu2]==1) ? -13 : 13;
      lept1_d0 = muD0[indexMu2];
      lept1_dz = muDz[indexMu2];
      lept1_normalizedChi2 = muChi2NDF[indexMu2];
      lept1_nValidMuonHits = muMuonHits[indexMu2];
      lept1_nMatchedStations = muStations[indexMu2];
      lept1_nValidPixelHits = muPixelHits[indexMu2];
      lept1_nTrackerLayers = muTrkLayers[indexMu2];
      lept1_muonBestTrack_dxyVTX = muD0[indexMu2];
      lept1_muonBestTrack_dzVTX = muDz[indexMu2];
      lept1_ptError = muBestTrkPtError[indexMu2];
      lept1_chiso = muPFChIso03[indexMu2];
      lept1_phoiso = muPFPhoIso03[indexMu2];
      lept1_neuiso = muPFNeuIso03[indexMu2];
      lept1_SIP = muSIP[indexMu2];
      lept1_sigEOverE = muBestTrkPtError[indexMu2]/muBestTrkPt[indexMu2];
      */
      if (mc) {
	lept0_SelSF = muSF(lept0_pt, lept0_eta);
	lept0_trigSF = muTrigSF_leg1(lept0_pt, lept0_eta);
	lept0_SelSF_Err = muSF_err(lept0_pt, lept0_eta);
	lept0_trigSF_Err = muTrigSF_leg1_err(lept0_pt, lept0_eta);

	lept1_SelSF = muSF(lept1_pt, lept1_eta);
	lept1_trigSF = muTrigSF_leg2(lept1_pt, lept1_eta);
	lept1_SelSF_Err = muSF_err(lept1_pt, lept1_eta);
	lept1_trigSF_Err = muTrigSF_leg2_err(lept1_pt, lept1_eta);
	lept_dzSF = mu_trigSF_DZ(abs(lept0_eta), abs(lept1_eta));
	lept_dzSF_Err = mu_trigSF_DZ_err(abs(lept0_eta), abs(lept1_eta));

	//cout << "Err of Sel SF for mu1: " << lept0_SelSF_Err << "\t mu2: " << lept1_SelSF_Err << endl;
	if ( lept0_SelSF_Err <0.000001) cout << lept0_pt << "\t eta: " << lept0_eta << endl;


      }

      //cout << "z info from muon channel" << endl;
      pair = mu[0] + mu[1];
      z_mass = pair.M();
      z_pt = pair.Pt();
      z_y = pair.Rapidity();
      z_eta = pair.Eta();
      z_phi = pair.Phi();
      z_charge = muCh[indexMu1] + muCh[indexMu2];

      deltaPhi_lept = deltaPhi(mu[0].Phi(),mu[1].Phi());
      deltaR_lept = deltaR(mu[0].Eta(), mu[0].Phi(),mu[1].Eta(), mu[1].Phi());

      if (pair.M() > 60 && pair.M() < 120) {
	if (z_charge == 0) heemass->Fill(pair.M());
	else heemassSS->Fill(pair.M());
      } 
    }
    else {
      z_mass = -99.;
      z_pt  = -99.;
      z_y = -99.;
      z_eta = -99.;
      z_phi = -99.;
      z_charge = -99;
    }
    
    //cout << "z_mass for WWToLNuQQ: " << z_mass << endl;
    if (z_mass < 0.) continue;

    pair_dPhi = -999.;
    //L1_EMTF problem of double muon trigger
    if (nmu > 1) {
      float pair_dPhi_ = deltaPhi(mu[0].Phi(), mu[1].Phi());
      if ( mu[0].Eta()*mu[1].Eta() > 0 && fabs(mu[0].Eta()) > 1.2 && fabs(mu[1].Eta()) > 1.2 ) 
	pair_dPhi = pair_dPhi_ * 180/TMath::Pi();
      else pair_dPhi = 999.;
    }
	
    //remove double counting for leptons in both dataset
    if (inputfile.Contains("DoubleEG") ) {
      if (leptType == 13 || (trig_Mu17_Mu8 == 1 && !trig_Ele23_Ele12)) continue;
    }
    if (inputfile.Contains("DoubleMu") ) {
      if (leptType == 11 || trig_Ele23_Ele12 == 1) continue;
    }

    //tZ->Fill();

    if ( nPho < 1) continue;
    npho = 0;

    vector<int> PreselectedPho;
    preselect_pho_Zg(PreselectedPho);

    if ( PreselectedPho.size() < 1) continue;

    for (vector<int>::const_iterator ipho = PreselectedPho.begin(); ipho != PreselectedPho.end(); ipho++) {
      //matching to gen level
      if ( mc ) {
	//if ( inputfile.Contains("Zg") || inputfile.Contains("ZGToLLG") || inputfile.Contains("TTGamma")) {
	if ( inputfile.Contains("Zg") || inputfile.Contains("ZGToLLG")) {
	  if ( !phoMatcher_prompt(data1, *ipho) ) continue;
	} else {
	  if (phoMatcher_prompt(data1, *ipho)) continue; //jet fake photon
	  //if (!phoMatch_ele(data1,*ipho)) continue; //electron fake photon
	  //if ( !phoMatch_jet(data1, *ipho) ) continue; //remove jet fake
	}
      }

      //phoCalibEt[*ipho] *= phoScale_rho_up[*ipho];
      //TRandom3 *gen = new TRandom3(0);
      //phoCalibEt[*ipho] *= (gen->Gaus(1,phoResol_rho_up[*ipho]));
      if (phoCalibEt[*ipho] < 2. ) continue;

      //cout << "check dR with lepton" << endl;
      float dR1, dR2;
      float dphi1 = -99, dphi2 = -99.;

      if (leptType == 11) {
	dR1 = deltaR(ele[0].Eta(),ele[0].Phi(),phoEta[*ipho],phoPhi[*ipho]);
	dR2 = deltaR(ele[1].Eta(),ele[1].Phi(),phoEta[*ipho],phoPhi[*ipho]);
	dphi1 = deltaPhi(ele[0].Phi(),phoPhi[*ipho]);
	dphi2 =	deltaPhi(ele[1].Phi(),phoPhi[*ipho]);
      }
      else if (leptType == 13) {
	dR1 = deltaR(mu[0].Eta(),mu[0].Phi(),phoEta[*ipho],phoPhi[*ipho]);
	dR2 = deltaR(mu[1].Eta(),mu[1].Phi(),phoEta[*ipho],phoPhi[*ipho]);
	dphi1 = deltaPhi(mu[0].Phi(),phoPhi[*ipho]);
	dphi2 =	deltaPhi(mu[1].Phi(),phoPhi[*ipho]);
      }
      else {
	dR1 = 99.;
	dR2 = 99.;
      }

      if (dR1 < 0.7 || dR2 < 0.7) continue;  //for SM_Zg

      TLorentzVector Photon, XZG;
      Photon.SetPtEtaPhiM(phoCalibEt[*ipho], phoEta[*ipho], phoPhi[*ipho], 0);
      XZG = pair + Photon;

      gamma_SF = 1.;
      gamma_SF_Err = 0.;
      gamma_CSEV_SF = 1.;
      gamma_CSEV_SF_Err = 0.;

      if ( mc ) {
	gamma_SF = phoSF( phoCalibEt[*ipho], phoSCEta[*ipho]);
	gamma_SF_Err = phoSF_err(phoCalibEt[*ipho], phoSCEta[*ipho]);
	gamma_CSEV_SF = csev_SF(phoCalibEt[*ipho], phoSCEta[*ipho]);
	gamma_CSEV_SF_Err = csev_SF_err(phoCalibEt[*ipho], phoSCEta[*ipho]);	
      }

      //cout << "fill photon variables" << endl;
      if (npho == 0) {
	if ( fabs(phoSCEta[*ipho]) < 1.4442) isEB = 1;
	else isEE = 1;
	gamma_pt = phoCalibEt[*ipho];
	gamma_eta = phoEta[*ipho];
	gamma_sceta = phoSCEta[*ipho];
	gamma_phi = phoPhi[*ipho];
	gamma_HoverE = phoHoverE[*ipho];
	gamma_mva = phoIDMVA[*ipho];
	gamma_ChWorstIso = phoPFChWorstIso[*ipho];
	gamma_PixSeed = phohasPixelSeed[*ipho];
	gamma_EleVeto = phoEleVeto[*ipho];
	//gamma_EnErr = phoEnErr[*ipho];
	gamma_En = phoE[*ipho];
	gamma_ssmva = PhotonSSMVA(data1, *ipho, tgr);
	gamma_ssmva_norw = PhotonSSMVA(data1, *ipho);
    
	/*
	etawei = 1;
	if (TMath::Abs(phoSCEta[*ipho])>1.5) {
	  if (phoSCEta[ipho]>-2.5 && phoSCEta[ipho]<-2.3)  etawei = 2.422;
	  else if ( phoSCEta[ipho]>-2.3 && phoSCEta[ipho]<-2.1) etawei = 1.194;
	  else if ( phoSCEta[ipho]>-2.1 && phoSCEta[ipho]<-1.9) etawei = 1.343;
	  else if ( phoSCEta[ipho]>-1.9 && phoSCEta[ipho]<-1.7) etawei = 0.617;
	  else if ( phoSCEta[ipho]>-1.7 && phoSCEta[ipho]<-1.5) etawei = 0.474;
	  else if ( phoSCEta[ipho]> 1.5 && phoSCEta[ipho]< 1.7) etawei = 0.678;
	  else if ( phoSCEta[ipho]> 1.7 && phoSCEta[ipho]< 1.9) etawei = 0.516;
	  else if ( phoSCEta[ipho]> 1.9 && phoSCEta[ipho]< 2.1) etawei = 1.070;
	  else if ( phoSCEta[ipho]> 2.1 && phoSCEta[ipho]< 2.3) etawei = 1.444;
	  else if ( phoSCEta[ipho]> 2.1 && phoSCEta[ipho]< 2.5) etawei = 2.825;
	  }*/
    
	gammaScale_stat_up = phoScale_stat_up[*ipho];
	gammaScale_stat_dn = phoScale_stat_dn[*ipho];
	gammaScale_syst_up = phoScale_syst_up[*ipho];
	gammaScale_syst_dn = phoScale_syst_dn[*ipho];
	gammaScale_gain_up = phoScale_gain_up[*ipho];
	gammaScale_gain_dn = phoScale_gain_dn[*ipho];
	gammaResol_rho_up = phoResol_rho_up[*ipho];
	gammaResol_rho_dn = phoResol_rho_dn[*ipho];
	gammaResol_phi_up = phoResol_phi_up[*ipho];
	gammaResol_phi_dn = phoResol_phi_dn[*ipho];

      

	//rho correction charged iso
	float chEA[] = {0.0360, 0.0377, 0.0306, 0.0283, 0.0254, 0.0217, 0.0167};
	float neuEA[] = {0.0597, 0.0807, 0.0629, 0.0197, 0.0184, 0.0284, 0.0591}; 
	float phoEA[] = {0.1210, 0.1107, 0.0699, 0.1056, 0.1457, 0.1719, 0.1998};
	int bin;
	if (phoSCEta[*ipho] < 1.0) bin = 0;
	else if (phoSCEta[*ipho] >= 1.0 && phoSCEta[*ipho] < 1.479) bin = 1;
	else if (phoSCEta[*ipho] >= 1.479 && phoSCEta[*ipho] < 2.0) bin = 2;
	else if (phoSCEta[*ipho] >= 2.0 && phoSCEta[*ipho] < 2.2) bin = 3;
	else if (phoSCEta[*ipho] >= 2.2 && phoSCEta[*ipho] < 2.3) bin = 4;
	else if (phoSCEta[*ipho] >= 2.3 && phoSCEta[*ipho] < 2.4) bin = 5;
	else bin = 6; // phoSCEta >= 2.4
	double corr = phoPFChIso[*ipho] - chEA[bin] * rho;
	double neucorr = phoPFNeuIso[*ipho] - neuEA[bin] * rho;
	double phocorr = phoPFPhoIso[*ipho] - phoEA[bin] * rho;
	
	//gamma_ChIso = phoPFChIso[*ipho];
	gamma_ChIso = TMath::Max(corr, 0.);
	gamma_NeuIso = TMath::Max(neucorr, 0.);
	gamma_PhoIso = TMath::Max(phocorr, 0.);

	//showershape correction
	gamma_SCEtaWidth = phoSCEtaWidth[*ipho];
	gamma_s4Full5x5 = phoE2x2Full5x5[*ipho]/phoE5x5Full5x5[*ipho];
	gamma_R9 = phoR9[*ipho];
	gamma_sigmaIetaIeta = phoSigmaIEtaIEtaFull5x5[*ipho];
	gamma_sigmaIetaIphi = phoSigmaIEtaIPhiFull5x5[*ipho];
	gamma_SCRawE = phoSCRawE[*ipho];
	gamma_scphiwidth = phoSCPhiWidth[*ipho];
	gamma_ESEnP1 = phoESEnP1[*ipho];
	gamma_ESEnP2 = phoESEnP2[*ipho];
	gamma_esRR = phoESEffSigmaRR[*ipho];

	if (data1.HasMC()) {
	  if(TMath::Abs(phoSCEta[*ipho])<1.5) {
	    gamma_SCEtaWidth_rw    = tgr[0]->Eval(gamma_SCEtaWidth);
	    gamma_s4Full5x5_rw     = tgr[1]->Eval(gamma_s4Full5x5);
	    gamma_R9_rw            = tgr[2]->Eval(gamma_R9);
	    gamma_sigmaIetaIeta_rw = tgr[3]->Eval(gamma_sigmaIetaIeta);
	    gamma_sigmaIetaIphi_rw = tgr[4]->Eval(gamma_sigmaIetaIphi);
	  } else {
	    gamma_SCEtaWidth_rw    = tgr[7]->Eval(gamma_SCEtaWidth);
	    gamma_s4Full5x5_rw     = tgr[8]->Eval(gamma_s4Full5x5);
	    gamma_R9_rw            = tgr[9]->Eval(gamma_R9);
	    gamma_sigmaIetaIeta_rw = tgr[10]->Eval(gamma_sigmaIetaIeta);
	    gamma_sigmaIetaIphi_rw = tgr[11]->Eval(gamma_sigmaIetaIphi);
	    
	  }
	}

	deltaR_PhoLept0 = dR1;
	deltaR_PhoLept1 = dR2;
	deltaPhi_PhoLept0 = dphi1;
	deltaPhi_PhoLept1 = dphi2;
	boss_mass = XZG.M();
	boss_pt  = XZG.Pt();
	boss_eta = XZG.Eta();
	boss_phi = XZG.Phi();
	//boss_massrel = sqrt( pow(lept0_EnErr/lept0_En, 2) + pow(lept1_EnErr/lept1_En, 2) + pow(gamma_EnErr/gamma_En,2) );
	//outtree->Fill();
      }
      npho++;
    }

    //loop on jet
    int njetmax = 40;
    //initial
    for (int i = 0; i < njetmax; i++) {
      jetPt_[i] = -99.;
      jetEta_[i] = -99.;
      dRLep1Jet[i] = -99.;
      dRLep2Jet[i] = -99.;
      dRPhoJet[i] = -99.;
    }


    //cout << "jet looop" << endl;
    //TLorentzVector selectedjet[njetmax];
    nselJet = 0;
    for (int ij = 0; ij < nJet; ij++) {
      if (data1.HasMC()) jetPt[ij] *= jetP4SmearUp[ij];
      //if (data1.HasMC()) jetPt[ij] *= jetP4SmearDo[ij];
      //jetPt[ij] *= (1-jetJECUnc[ij]);
      //cout << "check and eta" << endl;
      if (jetPt[ij] < 30. ) continue;
      if (abs(jetEta[ij]) > 2.4) continue;
      bool looseID = false;
      if ( abs(jetEta[ij]) <= 2.7) 
	looseID = (jetNHF[ij]<0.99 && jetNEF[ij]<0.99 && (jetNCH[ij]+jetNNP[ij])>1 )&& ((abs(jetEta[ij])<=2.4 && jetCHF[ij]>0 && jetNCH[ij]>0 && jetCEF[ij]<0.99) || abs(jetEta[ij])>2.4);
      else if (  abs(jetEta[ij]) <= 3.0)
	looseID = (jetNEF[ij]>0.01 && jetNHF[ij]<0.98 && jetNNP[ij]>2);
      else looseID = (jetNEF[ij]<0.90 && jetNNP[ij]>10);

      if (!looseID) continue;
	  
      //cout  <<" deltaR between jet and other" << endl;
      //separation between jet and lepton/photon    
      if (leptType == 11) {
	if (deltaR(jetEta[ij], jetPhi[ij], ele[0].Eta(),ele[0].Phi()) < 0.4) continue;
	if (deltaR(jetEta[ij], jetPhi[ij], ele[1].Eta(),ele[1].Phi()) < 0.4) continue;
      }
      else if (leptType == 13) {
	if (deltaR(jetEta[ij], jetPhi[ij], mu[0].Eta(),mu[0].Phi()) < 0.4) continue;
	if (deltaR(jetEta[ij], jetPhi[ij], mu[1].Eta(),mu[1].Phi()) < 0.4) continue;
      }

      if (gamma_pt > 0.) 
	if (deltaR(jetEta[ij], jetPhi[ij], gamma_eta, gamma_phi) < 0.4) continue;

      
      jetPt_[nselJet] = jetPt[ij];
      jetEta_[nselJet] = jetEta[ij];
      //jetP4Smear_[nselJet] = jetP4Smear[ij];
      //jetP4SmearUp_[nselJet] = jetP4SmearUp[ij];
      //jetP4SmearDo_[nselJet] = jetP4SmearDo[ij];
      //jetJECUnc_[nselJet] = jetJECUnc[ij];

      if (leptType == 11 ) {
	dRLep1Jet[nselJet] = deltaR(jetEta[ij], jetPhi[ij], ele[0].Eta(),ele[0].Phi());
	dRLep2Jet[nselJet] = deltaR(jetEta[ij], jetPhi[ij], ele[1].Eta(),ele[1].Phi());
      }
      else {
	dRLep1Jet[nselJet] = deltaR(jetEta[ij], jetPhi[ij], mu[0].Eta(),mu[0].Phi());
	dRLep2Jet[nselJet] = deltaR(jetEta[ij], jetPhi[ij], mu[1].Eta(),mu[1].Phi());
      }
      dRPhoJet[nselJet] = deltaR(jetEta[ij], jetPhi[ij], gamma_eta, gamma_phi);
      nselJet++;

    } //end of jet loop
    
    outtree->Fill();
    
  }

  hntotweight->SetBinContent(1, ntot);
  cout << "total of effective events: " << ntot << endl;
  cout << "total of photon with status flag = 0 is: " << nPhoPrompt << endl;
  cout << "total of number of gen having 2 leptons: " << ntotLep << endl;
  cout << "total positive events having 2 leptons + 1 photon: " << totpos << endl;
  cout << "total negative events having 2 leptons + 1 photon: " << totneg << endl;
  fo->cd();
  outtree->Write("", TObject::kOverwrite);
  fo->Write();
  fo->Close();
  delete fo;

  //frunlumi.close();

}
