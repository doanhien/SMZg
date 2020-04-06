#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

void makeMatrix_isPVGood(bool lhe = true, bool scale = 1, int num = 0, int iter = 0) {

  cout << "running " << iter << " th times "<< endl;

  TFile *fmc = new TFile("minitrees/ZGToLLG_5f_NoPtCut.root", "read");
  //TFile *fmc = new TFile("minitrees/ZGToLLG_LO_Madgraph_NoPtCut.root", "read");

  TTree *tmc = (TTree*) fmc->Get("outtree");

  if(!tmc) {
    cout<<"could not read MC tree\n";
  }


  float gamma_pt, gamma_eta, gamma_phi;
  float gamma_ssmva, gamma_ChIso;
  int isEB, isEE;
  float lept0_pt, lept0_eta, lept0_phi;
  float lept1_pt, lept1_eta, lept1_phi;
  int leptType;
  float pair_dPhi;
  int trig_Ele23_Ele12, trig_Mu17_Mu8;
  float z_mass, boss_mass;
  float boss_pt;
  int z_charge;
  float puweigj_65nb, genWeight;
  float lept0_RecoSF, lept1_RecoSF;
  float lept0_SelSF, lept1_SelSF;
  float gamma_SF, gamma_CSEV_SF;
  float lept0_trigSF, lept1_trigSF;
  float lept_dzSF;
  int nselJet;
  float jetEta[10], jetPt[10];
  float jetJECUnc[10];
  float jetP4Smear[10], jetP4SmearUp[10], jetP4SmearDo[10];

  float gammaScale_stat_up;
  float gammaScale_stat_dn;
  float gammaScale_syst_up;
  float gammaScale_syst_dn;
  float gammaScale_gain_up;
  float gammaScale_gain_dn;
  float gammaResol_rho_up;
  float gammaResol_rho_dn;
  float gammaResol_phi_up;
  float gammaResol_phi_dn;

  float eleScale_stat_up_1;
  float eleScale_stat_dn_1;
  float eleScale_syst_up_1;
  float eleScale_syst_dn_1;
  float eleScale_gain_up_1;
  float eleScale_gain_dn_1;
  float eleResol_rho_up_1;
  float eleResol_rho_dn_1;
  float eleResol_phi_up_1;
  float eleResol_phi_dn_1;

  float eleScale_stat_up_2;
  float eleScale_stat_dn_2;
  float eleScale_syst_up_2;
  float eleScale_syst_dn_2;
  float eleScale_gain_up_2;
  float eleScale_gain_dn_2;
  float eleResol_rho_up_2;
  float eleResol_rho_dn_2;
  float eleResol_phi_up_2;
  float eleResol_phi_dn_2;


  float genPhoEt, genPhoEta;
  float genlepPt[2], genlepEta[2];
  float genZm, genXm;
  float genXpt;
  int ngenjet;


  //lhe level variables
  float lheEle1, lheEle2;
  float lheEleEta1, lheEleEta2;
  float lheMu1, lheMu2;
  float lheMuEta1, lheMuEta2;
  float lhePho1, lhePhoEta1;
  float lhe_Zm, lhe_Xm;
  float dRlhePhoEle1, dRlhePhoEle2;
  float dRlhePhoMu1, dRlhePhoMu2;
  int lhe_type;


  tmc->ResetBranchAddresses();
  tmc->SetBranchAddress("gamma_pt",         &gamma_pt);
  tmc->SetBranchAddress("gamma_eta",        &gamma_eta);
  tmc->SetBranchAddress("gamma_phi",        &gamma_phi);
  tmc->SetBranchAddress("gamma_ssmva",      &gamma_ssmva);
  tmc->SetBranchAddress("gamma_ChIso",      &gamma_ChIso);
  tmc->SetBranchAddress("isEB",             &isEB);
  tmc->SetBranchAddress("isEE",             &isEE);
  tmc->SetBranchAddress("lept0_pt",         &lept0_pt);
  tmc->SetBranchAddress("lept0_eta",        &lept0_eta);
  tmc->SetBranchAddress("lept0_phi",        &lept0_phi);
  tmc->SetBranchAddress("lept1_pt",         &lept1_pt);
  tmc->SetBranchAddress("lept1_eta",        &lept1_eta);
  tmc->SetBranchAddress("lept1_phi",        &lept1_phi);
  tmc->SetBranchAddress("leptType",         &leptType);
  tmc->SetBranchAddress("trig_Ele23_Ele12", &trig_Ele23_Ele12);
  tmc->SetBranchAddress("trig_Mu17_Mu8",    &trig_Mu17_Mu8);
  tmc->SetBranchAddress("z_mass",           &z_mass);
  tmc->SetBranchAddress("boss_mass",        &boss_mass);
  tmc->SetBranchAddress("boss_pt",          &boss_pt);
  tmc->SetBranchAddress("z_charge",         &z_charge);
  tmc->SetBranchAddress("puweigj_65nb",     &puweigj_65nb);
  tmc->SetBranchAddress("genWeight",        &genWeight);
  tmc->SetBranchAddress("lept0_RecoSF",     &lept0_RecoSF);
  tmc->SetBranchAddress("lept1_RecoSF",     &lept1_RecoSF);
  tmc->SetBranchAddress("lept0_SelSF",      &lept0_SelSF);
  tmc->SetBranchAddress("lept1_SelSF",      &lept1_SelSF);
  tmc->SetBranchAddress("gamma_SF",         &gamma_SF);
  tmc->SetBranchAddress("gamma_CSEV_SF",    &gamma_CSEV_SF);
  tmc->SetBranchAddress("lept0_trigSF",     &lept0_trigSF);
  tmc->SetBranchAddress("lept1_trigSF",     &lept1_trigSF);
  tmc->SetBranchAddress("lept_dzSF",        &lept_dzSF);
  tmc->SetBranchAddress("nselJet",          &nselJet);
  tmc->SetBranchAddress("jetEta",           &jetEta);
  tmc->SetBranchAddress("jetPt",            &jetPt);
  tmc->SetBranchAddress("jetP4Smear",       &jetP4Smear);
  tmc->SetBranchAddress("jetP4SmearUp",     &jetP4SmearUp);
  tmc->SetBranchAddress("jetP4SmearDo",     &jetP4SmearDo);
  tmc->SetBranchAddress("jetJECUnc",        &jetJECUnc);
  tmc->SetBranchAddress("pair_dPhi",        &pair_dPhi);

  tmc->SetBranchAddress("gammaScale_stat_up", &gammaScale_stat_up);
  tmc->SetBranchAddress("gammaScale_stat_dn", &gammaScale_stat_dn);
  tmc->SetBranchAddress("gammaScale_syst_up", &gammaScale_syst_up);
  tmc->SetBranchAddress("gammaScale_syst_dn", &gammaScale_syst_dn);
  tmc->SetBranchAddress("gammaScale_gain_up", &gammaScale_gain_up);
  tmc->SetBranchAddress("gammaScale_gain_dn", &gammaScale_gain_dn);
  tmc->SetBranchAddress("gammaResol_rho_up",  &gammaResol_rho_up);
  tmc->SetBranchAddress("gammaResol_rho_dn",  &gammaResol_rho_dn);
  tmc->SetBranchAddress("gammaResol_phi_up",  &gammaResol_phi_up);
  tmc->SetBranchAddress("gammaResol_phi_dn",  &gammaResol_phi_dn);

  tmc->SetBranchAddress("eleScale_stat_up_1", &eleScale_stat_up_1);
  tmc->SetBranchAddress("eleScale_stat_dn_1", &eleScale_stat_dn_1);
  tmc->SetBranchAddress("eleScale_syst_up_1", &eleScale_syst_up_1);
  tmc->SetBranchAddress("eleScale_syst_dn_1", &eleScale_syst_dn_1);
  tmc->SetBranchAddress("eleScale_gain_up_1", &eleScale_gain_up_1);
  tmc->SetBranchAddress("eleScale_gain_dn_1", &eleScale_gain_dn_1);
  tmc->SetBranchAddress("eleResol_rho_up_1",  &eleResol_rho_up_1);
  tmc->SetBranchAddress("eleResol_rho_dn_1",  &eleResol_rho_dn_1);
  tmc->SetBranchAddress("eleResol_phi_up_1",  &eleResol_phi_up_1);
  tmc->SetBranchAddress("eleResol_phi_dn_1",  &eleResol_phi_dn_1);

  tmc->SetBranchAddress("eleScale_stat_up_2", &eleScale_stat_up_2);
  tmc->SetBranchAddress("eleScale_stat_dn_2", &eleScale_stat_dn_2);
  tmc->SetBranchAddress("eleScale_syst_up_2", &eleScale_syst_up_2);
  tmc->SetBranchAddress("eleScale_syst_dn_2", &eleScale_syst_dn_2);
  tmc->SetBranchAddress("eleScale_gain_up_2", &eleScale_gain_up_2);
  tmc->SetBranchAddress("eleScale_gain_dn_2", &eleScale_gain_dn_2);
  tmc->SetBranchAddress("eleResol_rho_up_2",  &eleResol_rho_up_2);
  tmc->SetBranchAddress("eleResol_rho_dn_2",  &eleResol_rho_dn_2);
  tmc->SetBranchAddress("eleResol_phi_up_2",  &eleResol_phi_up_2);
  tmc->SetBranchAddress("eleResol_phi_dn_2",  &eleResol_phi_dn_2);

  tmc->SetBranchAddress("genPhoEt",         &genPhoEt);
  tmc->SetBranchAddress("genPhoEta",        &genPhoEta);
  tmc->SetBranchAddress("genlepPt",         genlepPt);
  tmc->SetBranchAddress("genlepEta",        genlepEta);
  tmc->SetBranchAddress("genZm",            &genZm);
  tmc->SetBranchAddress("genXm",            &genXm);
  tmc->SetBranchAddress("genXpt",           &genXpt);
  tmc->SetBranchAddress("ngenjet",          &ngenjet);

  tmc->SetBranchAddress("lheEle1",          &lheEle1);
  tmc->SetBranchAddress("lheEle2",          &lheEle2);
  tmc->SetBranchAddress("lheEleEta1",       &lheEleEta1);
  tmc->SetBranchAddress("lheEleEta2",       &lheEleEta2);
  tmc->SetBranchAddress("lheMu1",           &lheMu1);
  tmc->SetBranchAddress("lheMu2",           &lheMu2);
  tmc->SetBranchAddress("lheMuEta1",        &lheMuEta1);
  tmc->SetBranchAddress("lheMuEta2",        &lheMuEta2);
  tmc->SetBranchAddress("lhePho1",          &lhePho1);
  tmc->SetBranchAddress("lhePhoEta1",       &lhePhoEta1);
  tmc->SetBranchAddress("lhe_Zm",           &lhe_Zm);
  tmc->SetBranchAddress("lhe_Xm",           &lhe_Xm);
  tmc->SetBranchAddress("lhe_type",         &lhe_type);
  tmc->SetBranchAddress("dRlhePhoEle1",     &dRlhePhoEle1);
  tmc->SetBranchAddress("dRlhePhoEle2",     &dRlhePhoEle2);
  tmc->SetBranchAddress("dRlhePhoMu1",      &dRlhePhoMu1);
  tmc->SetBranchAddress("dRlhePhoMu2",      &dRlhePhoMu2);


  //binnings and histograms
  const int bin_gen = 11;
  //float pt_gen[bin_gen+1] = {15, 20, 25, 30, 35, 45, 55, 65, 75, 85, 95, 120, 1000};
  float pt_gen[bin_gen+1] = {20, 25, 30, 35, 45, 55, 65, 75, 85, 95, 120, 1000};

  //const int bin_rec = 25;
  //float pt_rec[bin_rec+1] = {15, 17, 19, 21, 23, 25, 27, 29, 32, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 105, 120, 200, 1000};
  const int bin_rec = 22;
  float pt_rec[bin_rec+1] = {20, 22, 25, 27, 29, 32, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 105, 120, 200, 1000};

  //for response matrix
  TH2D *hrec_gen_ele_EB = new TH2D("hrec_gen_ele_EB", "hrec_gen_ele_EB", bin_gen, pt_gen, bin_rec, pt_rec);
  TH2D *hrec_gen_ele_EE = new TH2D("hrec_gen_ele_EE", "hrec_gen_ele_EE", bin_gen, pt_gen, bin_rec, pt_rec);

  TH2D *hrec_gen_mu_EB = new TH2D("hrec_gen_mu_EB", "hrec_gen_mu_EB", bin_gen, pt_gen, bin_rec, pt_rec);
  TH2D *hrec_gen_mu_EE = new TH2D("hrec_gen_mu_EE", "hrec_gen_mu_EE", bin_gen, pt_gen, bin_rec, pt_rec);


  //------------------ mass -------------------------//
  const int bin_gen_mllg = 11;
  float mllg_gen[bin_gen_mllg+1] = {70, 88, 95, 110, 135, 170, 210, 270, 350, 470, 640, 3000};

  const int bin_rec_mllg = 22;
  float mllg_rec[bin_rec_mllg+1] = {70, 85, 88, 92, 95, 100, 105, 110, 115, 120, 135, 150, 170, 190, 210, 240, 270, 300, 350, 400, 470, 640, 3000};

  TH2D *hrec_gen_mllg_EB_ele = new TH2D("hrec_gen_mllg_EB_ele", "hrec_gen_mllg_EB_ele", bin_gen_mllg, mllg_gen, bin_rec_mllg, mllg_rec);
  TH2D *hrec_gen_mllg_EE_ele = new TH2D("hrec_gen_mllg_EE_ele", "hrec_gen_mllg_EE_ele", bin_gen_mllg, mllg_gen, bin_rec_mllg, mllg_rec);

  TH2D *hrec_gen_mllg_EB_mu = new TH2D("hrec_gen_mllg_EB_mu", "hrec_gen_mllg_EB_mu", bin_gen_mllg, mllg_gen, bin_rec_mllg, mllg_rec);
  TH2D *hrec_gen_mllg_EE_mu = new TH2D("hrec_gen_mllg_EE_mu", "hrec_gen_mllg_EE_mu", bin_gen_mllg, mllg_gen, bin_rec_mllg, mllg_rec);


  //-----------------pT of Zg ----------------//
  const int bin_gen_ptllg = 8;
  float ptllg_gen[bin_gen_ptllg+1] = {0, 10, 20, 30, 40, 50, 70, 120, 1000};

  const int bin_rec_ptllg = 16;
  float ptllg_rec[bin_rec_ptllg+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 90, 120, 150, 1000};

  TH2D *hrec_gen_ptllg_EB_ele = new TH2D("hrec_gen_ptllg_EB_ele", "hrec_gen_ptllg_EB_ele", bin_gen_ptllg, ptllg_gen, bin_rec_ptllg, ptllg_rec);
  TH2D *hrec_gen_ptllg_EE_ele = new TH2D("hrec_gen_ptllg_EE_ele", "hrec_gen_ptllg_EE_ele", bin_gen_ptllg, ptllg_gen, bin_rec_ptllg, ptllg_rec);

  TH2D *hrec_gen_ptllg_EB_mu = new TH2D("hrec_gen_ptllg_EB_mu", "hrec_gen_ptllg_EB_mu", bin_gen_ptllg, ptllg_gen, bin_rec_ptllg, ptllg_rec);
  TH2D *hrec_gen_ptllg_EE_mu = new TH2D("hrec_gen_ptllg_EE_mu", "hrec_gen_ptllg_EE_mu", bin_gen_ptllg, ptllg_gen, bin_rec_ptllg, ptllg_rec);


  //-------- binning for jet multiplicity --------------//
  const int njet_recbin = 7;
  const int njet_genbin = 4;

  TH2D *hrec_gen_ele_njet = new TH2D("hrec_gen_ele_njet", "hrec_gen_ele_njet", njet_genbin, 0, njet_genbin, njet_recbin, 0, njet_recbin) ;
  TH2D *hrec_gen_mu_njet = new TH2D("hrec_gen_mu_njet", "hrec_gen_mu_njet", njet_genbin, 0, njet_genbin, njet_recbin, 0, njet_recbin);

  cout << "loop over mc tree" << endl;
  cout << "total events from tree: " << tmc->GetEntriesFast() << endl;

  int countEle = 0;
  int countMu = 0;
  int passEv = 0;

  float ev_njet[4][7] = {0};
  float ev_njet_mu[4][7] = {0};

  //Long64_t totEv = tmc->GetEntriesFast()*0.8;
  Long64_t totEv = tmc->GetEntriesFast();
  //for trainning
  for (Long64_t iev = 0; iev < totEv; iev++) {
    tmc->GetEntry(iev);
    if(tmc->GetEntry(iev)<=0) break;
    //if ( (iev % 2) != 0 ) continue;
    if (iev%20000 == 0) cout << "processing event: " << iev+1 << " th" << endl;
    /*
    if (scale) {
      if (num == 0) gamma_pt *= gammaScale_stat_up;
      if (num == 1) gamma_pt *= gammaScale_stat_dn;
      if (num == 2) gamma_pt *= gammaScale_syst_up;
      if (num == 3) gamma_pt *= gammaScale_syst_dn;
      if (num == 4) gamma_pt *= gammaScale_gain_up;
      if (num == 5) gamma_pt *= gammaScale_gain_dn;
      if (num == -1) gamma_pt *= 1;
    }
    else {
      TRandom *gen = new TRandom3(0);
      if (num == 0) gamma_pt *= (gen->Gaus(1,gammaResol_rho_up));
      if (num == 1) gamma_pt *= (gen->Gaus(1,gammaResol_rho_dn));
      if (num == 2) gamma_pt *= (gen->Gaus(1,gammaResol_phi_up));
      if (num == 3) gamma_pt *= (gen->Gaus(1,gammaResol_phi_dn));
      if (num == -1) gamma_pt *= 1;
    }
    */
        
    if (trig_Ele23_Ele12 ==1 && leptType==11) { //scale and smear for ele
      if (scale) {
	if (num == 0) {lept0_pt *= eleScale_stat_up_1; lept1_pt *= eleScale_stat_up_2;}
	if (num == 1) {lept0_pt *= eleScale_stat_dn_1; lept1_pt *= eleScale_stat_dn_2;}
	if (num == 2) {lept0_pt *= eleScale_syst_up_1; lept1_pt *= eleScale_syst_up_2;}
	if (num == 3) {lept0_pt *= eleScale_syst_dn_1; lept1_pt *= eleScale_syst_dn_2;}
	if (num == 4) {lept0_pt *= eleScale_gain_up_1; lept1_pt *= eleScale_gain_up_2;}
	if (num == 5) {lept0_pt *= eleScale_gain_dn_1; lept1_pt *= eleScale_gain_dn_2;}
	if (num == -1) {lept0_pt *= 1; lept1_pt *= 1;}
      }
      else {
	TRandom *gen1 = new TRandom3(0);
	//TRandom *gen2 = new TRandom3(0);
	if (num == 0) {lept0_pt *= (gen1->Gaus(1,eleResol_rho_up_1)); lept1_pt *= (gen1->Gaus(1,eleResol_rho_up_2));}
	if (num == 1) {lept0_pt *= (gen1->Gaus(1,eleResol_rho_dn_1)); lept1_pt *= (gen1->Gaus(1,eleResol_rho_dn_2));}
	if (num == 2) {lept0_pt *= (gen1->Gaus(1,eleResol_phi_up_1)); lept1_pt *= (gen1->Gaus(1,eleResol_phi_up_2));}
	if (num == 3) {lept0_pt *= (gen1->Gaus(1,eleResol_phi_dn_1)); lept1_pt *= (gen1->Gaus(1,eleResol_phi_dn_2));}
	if (num == -1) {lept0_pt *= 1; lept1_pt *= 1; }
      }
    }
        

    if (lhe) {
      if (lhePho1 < 15. || lhe_Xm < 0.) continue; //for lhe level
      if (lhe_Zm < 50.) continue;
    }
    else {
      if (genPhoEt < 0. || genXm < 0.) continue; //for gen level
      if ( fabs(genPhoEta)>2.5 ) continue;
      if ( fabs(genlepEta[0]) > 2.4 || fabs(genlepEta[1]) > 2.4 ) continue;
      if (genlepPt[0] < 25 || genlepPt[1] < 20) continue;
      if (genZm < 50.) continue;
    }

    if (z_charge != 0) continue;
    if (gamma_pt < 0.) continue;

    passEv++;
    if (gamma_pt < 15.) { gamma_SF = 1.; gamma_CSEV_SF = 1.;}

    float truth_Pt = genPhoEt;
    if (lhe) truth_Pt = lhePho1;

    float truth_mllg = genXm;
    if (lhe) truth_mllg = lhe_Xm;

    //electron
    if ( trig_Ele23_Ele12 ==1 && leptType ==11 ) {
      if (lhe) {
	if (dRlhePhoEle1 < 0.7 || dRlhePhoEle2 < 0.7) continue;
	if (fabs(lheEleEta1) > 2.4 || fabs(lheEleEta2) > 2.4) continue;
      }
      
      TLorentzVector ele1, ele2, pho, X;
      ele1.SetPtEtaPhiM(lept0_pt, lept0_eta, lept0_phi, 0.5*0.001);
      ele2.SetPtEtaPhiM(lept1_pt, lept1_eta, lept1_phi, 0.5*0.001);
      pho.SetPtEtaPhiM(gamma_pt, gamma_eta, gamma_phi, 0.);
      X = ele1 + ele2 + pho;
      float Xm = X.M();
      float Xpt = X.Pt();
      float Zm = (ele1 + ele2).M();
      if (Zm < 50.) continue;
      if (ele1.Pt() < 25.) continue;
      if (ele2.Pt() < 20.) continue;
      if (fabs(ele1.Eta())> 2.4) continue;
      if (fabs(ele2.Eta())> 2.4) continue;
      //cout << "boss_mass: " << boss_mass << "\t Xm: " << Xm << endl;

      float wei = puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF;
      float weigen = puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF*lept0_trigSF*lept1_trigSF*lept_dzSF;

      if (isEB && gamma_ChIso<2.) {
        if (gamma_pt >= 20.)
          hrec_gen_ele_EB->Fill(truth_Pt,gamma_pt,wei);
        else if (truth_Pt >= 20.) hrec_gen_ele_EB->Fill(truth_Pt, -1, weigen);
        if (gamma_pt >= 20.) {
          //if (Xm >= 70. && truth_mllg>=70.) hrec_gen_mllg_EB_ele->Fill(truth_mllg, Xm, wei);
          if (Xm >= 70.) hrec_gen_mllg_EB_ele->Fill(truth_mllg, Xm, wei);
          else if (truth_mllg >= 70. && Xm<70.)  hrec_gen_mllg_EB_ele->Fill(truth_mllg, -1, weigen);
	  if (Xpt >=0. && genXpt >=0. ) hrec_gen_ptllg_EB_ele->Fill(genXpt, Xpt, wei);
	  else if (Xpt <0. && genXpt >=0. ) hrec_gen_ptllg_EB_ele->Fill(genXpt, -1, wei);
        }
      }
      else if( isEE && gamma_ChIso<1.5)  {
        if (gamma_pt >= 20.)
          hrec_gen_ele_EE->Fill(truth_Pt, gamma_pt, wei);
	else if (truth_Pt >= 20.) hrec_gen_ele_EE->Fill(truth_Pt, -1, weigen);
        if ( gamma_pt >= 20.) {
          //if (Xm >= 70. && truth_mllg>=70.) hrec_gen_mllg_EE_ele->Fill(truth_mllg, Xm, wei);
          if (Xm >= 70.) hrec_gen_mllg_EE_ele->Fill(truth_mllg, Xm, wei);
          else if(truth_mllg >= 70. && Xm<70.)  hrec_gen_mllg_EE_ele->Fill(truth_mllg, -1, weigen);
	  if (Xpt >=0. && genXpt >=0. ) hrec_gen_ptllg_EE_ele->Fill(genXpt, Xpt, wei);
	  else if (Xpt <0. && genXpt >=0. ) hrec_gen_ptllg_EE_ele->Fill(genXpt, -1, wei);

        }
      }

      //for jet multiplicity
      if ( gamma_pt >= 20. && ((isEB && gamma_ChIso<2.) || (isEE && gamma_ChIso<1.5)) ) {
	if (nselJet == 0) {
	  if (ngenjet==0) ev_njet[0][0] += wei;
	  if (ngenjet==1) ev_njet[1][0] += wei;
	  if (ngenjet==2) ev_njet[2][0] += wei;
	  if (ngenjet>2) ev_njet[3][0] += wei;
	}
	else if (nselJet==1) {
	  if ( jetPt[0] > 30. && jetPt[0] < 50.) {
	    if (ngenjet==0) ev_njet[0][1] += wei;
	    if (ngenjet==1) ev_njet[1][1] += wei;
	    if (ngenjet==2) ev_njet[2][1] += wei;
	    if (ngenjet>2) ev_njet[3][1] += wei;
	  }
	  else if (jetPt[0] >= 50.) {
	    if (ngenjet==0) ev_njet[0][2] += wei;
	    if (ngenjet==1) ev_njet[1][2] += wei;
	    if (ngenjet==2) ev_njet[2][2] += wei;
	    if (ngenjet>2) ev_njet[3][2] += wei;
	  }
	}
	else if (nselJet==2) {
	  if ( jetPt[1] > 30. && jetPt[1] <50) {
	    if (ngenjet==0) ev_njet[0][3] += wei;
	    if (ngenjet==1) ev_njet[1][3] += wei;
	    if (ngenjet==2) ev_njet[2][3] += wei;
	    if (ngenjet>2) ev_njet[3][3] += wei;
	  }
	  else if (jetPt[1] >= 50.) {
	    if (ngenjet==0) ev_njet[0][4] += wei;
	    if (ngenjet==1) ev_njet[1][4] += wei;
	    if (ngenjet==2) ev_njet[2][4] += wei;
	    if (ngenjet>2) ev_njet[3][4] += wei;
	  }
	}
	else {
	  if ( jetPt[2] > 30. && jetPt[2] < 50.) {
	    if (ngenjet==0) ev_njet[0][5] += wei;
	    if (ngenjet==1) ev_njet[1][5] += wei;
	    if (ngenjet==2) ev_njet[2][5] += wei;
	    if (ngenjet>2) ev_njet[3][5] += wei;
	  }
	  else if (jetPt[2] >= 50.) {
	    if (ngenjet==0) ev_njet[0][6] += wei;
	    if (ngenjet==1) ev_njet[1][6] += wei;
	    if (ngenjet==2) ev_njet[2][6] += wei;
	    if (ngenjet>2) ev_njet[3][6] += wei;
	  }
	}
      }
    }
  
    //muon
    if (trig_Mu17_Mu8==1 && leptType==13 && pair_dPhi>70.) {
      if (lhe) {
        if (dRlhePhoMu1 < 0.7 || dRlhePhoMu2 < 0.7) continue;
        if (fabs(lheMuEta1) > 2.4 || fabs(lheMuEta2) > 2.4 ) continue;
      }

      TLorentzVector mu1, mu2, X, pho;
      mu1.SetPtEtaPhiM(lept0_pt, lept0_eta, lept0_phi, 105*0.001);
      mu2.SetPtEtaPhiM(lept1_pt, lept1_eta, lept1_phi, 105*0.001);
      pho.SetPtEtaPhiM(gamma_pt, gamma_eta, gamma_phi, 0.);

      X = mu1 + mu2 + pho;
      float Xm = X.M();
      float Xpt = X.Pt();
      if (truth_mllg > 3000.) cout << "truth mass: " << truth_mllg << "\t reco: " << Xm << "\t pho pt: " << gamma_pt << endl;

      if ( (mu1+mu2).M() < 50.) continue;
      if (mu1.Pt() < 25.) continue;
      if (mu2.Pt() < 20.) continue;

      float wei = puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;
      float weigen = puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;

      if (isEB && gamma_ChIso<2.) {
        if (gamma_pt >= 20.)
          hrec_gen_mu_EB->Fill(truth_Pt, gamma_pt, wei);
        else if (truth_Pt >= 20.)
          hrec_gen_mu_EB->Fill(truth_Pt, -1, weigen);

        if (gamma_pt >= 20.) {
          //if (Xm >= 70. && truth_mllg>=70.) hrec_gen_mllg_EB_mu->Fill(truth_mllg, Xm, wei);
          if (Xm >= 70.) hrec_gen_mllg_EB_mu->Fill(truth_mllg, Xm, wei);
          else if (truth_mllg >= 70. && Xm<70.)  hrec_gen_mllg_EB_mu->Fill(truth_mllg, -1, weigen);
	  if (Xpt >=0. && genXpt >=0. ) hrec_gen_ptllg_EB_mu->Fill(genXpt, Xpt, wei);
          else if (Xpt <0. && genXpt >=0. ) hrec_gen_ptllg_EB_mu->Fill(genXpt, -1, wei);
        }

      }
      else if( isEE && gamma_ChIso<1.5)  {
        if (gamma_pt >= 20.)
          hrec_gen_mu_EE->Fill(truth_Pt, gamma_pt, wei);
        else if (truth_Pt >= 20.)
          hrec_gen_mu_EE->Fill(truth_Pt, -1, weigen);

        if (gamma_pt >= 20.) {
          //if (Xm >= 70. && truth_mllg>=70.) hrec_gen_mllg_EE_mu->Fill(truth_mllg, Xm, wei);
          if (Xm >= 70.) hrec_gen_mllg_EE_mu->Fill(truth_mllg, Xm, wei);
          else if (truth_mllg >= 70. && Xm<70.)  hrec_gen_mllg_EE_mu->Fill(truth_mllg, -1, weigen);
	  if (Xpt >=0. && genXpt >=0. ) hrec_gen_ptllg_EE_mu->Fill(genXpt, Xpt, wei);
          else if (Xpt <0. && genXpt >=0. ) hrec_gen_ptllg_EE_mu->Fill(genXpt, -1, wei);
        }
      }

      //for jet multiplicity
      if ( gamma_pt >= 20. && ((isEB && gamma_ChIso<2.) || (isEE && gamma_ChIso<1.5)) ) {
	if (nselJet==0) {
	  if (ngenjet==0) ev_njet_mu[0][0] += wei;
	  if (ngenjet==1) ev_njet_mu[1][0] += wei;
	  if (ngenjet==2) ev_njet_mu[2][0] += wei;
	  if (ngenjet>2) ev_njet_mu[3][0] += wei;
	}
	else if (nselJet==1) {
	  if ( jetPt[0] > 30. && jetPt[0] < 50.) {
	    if (ngenjet==0) ev_njet_mu[0][1] += wei;
	    if (ngenjet==1) ev_njet_mu[1][1] += wei;
	    if (ngenjet==2) ev_njet_mu[2][1] += wei;
	    if (ngenjet>2) ev_njet_mu[3][1] += wei;
	  }
	  if (jetPt[0] >= 50.) {
	    if (ngenjet==0) ev_njet_mu[0][2] += wei;
	    if (ngenjet==1) ev_njet_mu[1][2] += wei;
	    if (ngenjet==2) ev_njet_mu[2][2] += wei;
	    if (ngenjet>2) ev_njet_mu[3][2] += wei;
	  }
	}
	else if (nselJet==2) {
	  if ( jetPt[1] > 30. && jetPt[1] < 50.) {
	    if (ngenjet==0) ev_njet_mu[0][3] += wei;
	    if (ngenjet==1) ev_njet_mu[1][3] += wei;
	    if (ngenjet==2) ev_njet_mu[2][3] += wei;
	    if (ngenjet>2) ev_njet_mu[3][3] += wei;
	  }
	  if (jetPt[1] >= 50.) {
	    if (ngenjet==0) ev_njet_mu[0][4] += wei;
	    if (ngenjet==1) ev_njet_mu[1][4] += wei;
	    if (ngenjet==2) ev_njet_mu[2][4] += wei;
	    if (ngenjet>2) ev_njet_mu[3][4] += wei;
	  }
	}
	else {
	  if ( jetPt[2] > 30. && jetPt[2] < 50.) {
	    if (ngenjet==0) ev_njet_mu[0][5] += wei;
	    if (ngenjet==1) ev_njet_mu[1][5] += wei;
	    if (ngenjet==2) ev_njet_mu[2][5] += wei;
	    if (ngenjet>2) ev_njet_mu[3][5] += wei;
	  }
	  if (jetPt[2] >= 50.) {
	    if (ngenjet==0) ev_njet_mu[0][6] += wei;
	    if (ngenjet==1) ev_njet_mu[1][6] += wei;
	    if (ngenjet==2) ev_njet_mu[2][6] += wei;
	    if (ngenjet>2) ev_njet_mu[3][6] += wei;
	  }
	}
      }
    }
  }

  cout << "total event passing: " << passEv << endl;
  //cout << "testing histogram" << endl;

  //filling response matrix for jet multiplicity
  for (int i = 0; i < njet_genbin; i++) {
    for (int j = 0; j < njet_recbin; j++) {
      hrec_gen_ele_njet->SetBinContent(i+1, j+1 , ev_njet[i][j]);
      hrec_gen_mu_njet->SetBinContent(i+1, j+1 , ev_njet_mu[i][j]);
      cout << "ngen = " << i << "\t nrec=" << j << "\t" << ev_njet[i][j] << endl;
    }
  }


  //for output
  TString outname = "";
  if (lhe) outname = "output_EGM_Resol/ZgToLLg_5f_responseMatrix_isPVGood";
  else outname = "inputUnfold/SB_EE7to13_EE6to14/main/ZgToLLg_5f_responseMatrix_isPVGood_ele";
  //else outname = "inputUnfold/SB_EE7to13_EE6to14/EGM_Resol/ZgToLLg_5f_responseMatrix_isPVGood_ele";
  //else outname = "inputUnfold/SB_EE7to13_EE6to14/main/ZgToLLg_5f_responseMatrix_isPVGood_gencut";
  if (scale) {
    if (num == 0) outname+= "_scale_stat_up.root";
    else if (num == 1) outname+= "_scale_stat_dn.root";
    else if (num == 2) outname+= "_scale_syst_up.root";
    else if (num == 3) outname+= "_scale_syst_dn.root";
    else if (num == 4) outname+= "_scale_gain_up.root";
    else if (num == 5) outname+= "_scale_gain_dn.root";
    else if (num == -1) outname+= ".root";
  }
  else {
    if (num == 0) outname+= "_resol_rho_up";
    else if (num == 1) outname+= "_resol_rho_dn";
    else if (num == 2) outname+= "_resol_phi_up";
    else outname+= "_resol_phi_dn";
    outname += Form("_%d.root", iter);
  }


  TFile *fout = new TFile(outname, "recreate");
  fout->cd();

  //matrix
  hrec_gen_ele_EB->Write("hrec_gen_ele_EB");
  hrec_gen_ele_EE->Write("hrec_gen_ele_EE");
  hrec_gen_mu_EB->Write("hrec_gen_mu_EB");
  hrec_gen_mu_EE->Write("hrec_gen_mu_EE");
  hrec_gen_mllg_EB_ele->Write("hrec_gen_mllg_EB_ele");
  hrec_gen_mllg_EE_ele->Write("hrec_gen_mllg_EE_ele");
  hrec_gen_mllg_EB_mu->Write("hrec_gen_mllg_EB_mu");
  hrec_gen_mllg_EE_mu->Write("hrec_gen_mllg_EE_mu");
  hrec_gen_ele_njet->Write("hrec_gen_ele_njet");
  hrec_gen_mu_njet->Write("hrec_gen_mu_njet");

  hrec_gen_ptllg_EB_ele->Write("hrec_gen_ptllg_EB_ele");
  hrec_gen_ptllg_EE_ele->Write("hrec_gen_ptllg_EE_ele");
  hrec_gen_ptllg_EB_mu->Write("hrec_gen_ptllg_EB_mu");
  hrec_gen_ptllg_EE_mu->Write("hrec_gen_ptllg_EE_mu");

  fout->Write();
  fout->Close();
  fmc->Close();

}
