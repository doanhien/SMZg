//this is data yield from fitting
// on 25.12.18 with new sidbeband window
//EB: 7 < charged is < 13 GeV
//EE: 6 < charged iso < 14 GeV

#include "TTree.h"
#include "TH1.h"
#include "TH2.h"

void UnfoldInput(TString sample = "amcnlo") {

  TString fname;
  //if (sample.Contains("amcnlo")) fname = "../../ana_jet/fitting/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  if (sample.Contains("amcnlo")) fname = "minitrees/ZGToLLG_LO_Madgraph_NoPtCut.root";
  else if (sample.Contains("sherpa")) fname ="../../ana_jet/fitting/minitrees/ZGToLLG_5f_Summer16_TMVA420_UpTp6000_5VarCorr.root";
  TFile *fmc = new TFile(fname, "read");

  cout << "input file: " << fname << endl;

  TTree *tmc = (TTree*) fmc->Get("outtree");
  TH1F *htotwei = (TH1F*) fmc->Get("hntotweight");
  if(!tmc) {
    cout<<"could not read MC tree\n";
  }

  float gamma_pt, gamma_eta, gamma_phi;
  float gamma_ssmva, gamma_ChIso;
  int isEB, isEE;
  float lept0_pt, lept0_eta, lept0_phi;
  float lept1_pt, lept1_eta, lept1_phi;
  int leptType;
  int trig_Ele23_Ele12, trig_Mu17_Mu8;
  float z_mass, boss_mass;
  float boss_pt;
  int z_charge;
  float puweigj_65nb, genWeight;
  float lept0_RecoSF, lept1_RecoSF;
  float lept0_SelSF, lept1_SelSF;
  float gamma_SF;
  float gamma_CSEV_SF;
  int nselJet;
  float jetEta[10], jetPt[10];
  float pair_dPhi;
  float genPhoEt, genPhoEta;
  float genlepPt[2], genlepEta[2];
  float genZm, genXm, genXpt;
  int ngenjet;

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
  tmc->SetBranchAddress("nselJet",          &nselJet);
  tmc->SetBranchAddress("jetEta",           &jetEta);
  tmc->SetBranchAddress("jetPt",            &jetPt);
  tmc->SetBranchAddress("pair_dPhi",        &pair_dPhi);

  tmc->SetBranchAddress("genPhoEt",         &genPhoEt);
  tmc->SetBranchAddress("genPhoEta",        &genPhoEta);
  tmc->SetBranchAddress("genlepPt",         genlepPt);
  tmc->SetBranchAddress("genlepEta",        genlepEta);
  tmc->SetBranchAddress("genZm",            &genZm);
  tmc->SetBranchAddress("genXm",            &genXm);
  tmc->SetBranchAddress("genXpt",           &genXpt);
  tmc->SetBranchAddress("ngenjet",          &ngenjet);

  //binnings and histograms
  const int bin_gen = 11;
  float pt_gen[bin_gen+1] = {20, 25, 30, 35, 45, 55, 65, 75, 85, 95, 120, 1000};

  //inclusive
  const int bin_rec = 22;
  float pt_rec[bin_rec+1] = {20, 22, 25, 27, 29, 32, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 105, 120, 200, 1000};

  //for reco histogram
  TH1D *hrec_ele_EB = new TH1D("hrec_ele_EB", "hrec_ele_EB", bin_rec, pt_rec);
  TH1D *hrec_ele_EE = new TH1D("hrec_ele_EE", "hrec_ele_EE", bin_rec, pt_rec);

  TH1D *hrec_mu_EB = new TH1D("hrec_mu_EB", "hrec_mu_EB", bin_rec, pt_rec);
  TH1D *hrec_mu_EE = new TH1D("hrec_mu_EE", "hrec_mu_EE", bin_rec, pt_rec);

  //gen histogram for testing
  TH1D *hgen_ele_EB = new TH1D("hgen_ele_EB", "hgen_ele_EB", bin_gen, pt_gen);
  TH1D *hgen_ele_EE = new TH1D("hgen_ele_EE", "hgen_ele_EE", bin_gen, pt_gen);

  TH1D *hgen_mu_EB = new TH1D("hgen_mu_EB", "hgen_mu_EB", bin_gen, pt_gen);
  TH1D *hgen_mu_EE = new TH1D("hgen_mu_EE", "hgen_mu_EE", bin_gen, pt_gen);

  //------------------ mass -------------------------//

  const int bin_gen_mllg = 11;
  float mllg_gen[bin_gen_mllg+1] = {70, 88, 95, 110, 135, 170, 210, 270, 350, 470, 640, 3000};
  const int bin_rec_mllg = 22;
  float mllg_rec[bin_rec_mllg+1] = {70, 85, 88, 92, 95, 100, 105, 110, 115, 120, 135, 150, 170, 190, 210, 240, 270, 300, 350, 400, 470, 640, 3000};

  TH1D *hrec_mllg_EB_ele = new TH1D("hrec_mllg_EB_ele", "hrec_mllg_EB_ele", bin_rec_mllg, mllg_rec);
  TH1D *hrec_mllg_EE_ele = new TH1D("hrec_mllg_EE_ele", "hrec_mllg_EE_ele", bin_rec_mllg, mllg_rec);
  TH1D *hrec_mllg_EB_mu = new TH1D("hrec_mllg_EB_mu", "hrec_mllg_EB_mu", bin_rec_mllg, mllg_rec);
  TH1D *hrec_mllg_EE_mu = new TH1D("hrec_mllg_EE_mu", "hrec_mllg_EE_mu", bin_rec_mllg, mllg_rec);

  TH1D *hgen_mllg_EB_ele = new TH1D("hgen_mllg_EB_ele", "hgen_mllg_EB_ele", bin_gen_mllg, mllg_gen);
  TH1D *hgen_mllg_EE_ele = new TH1D("hgen_mllg_EE_ele", "hgen_mllg_EE_ele", bin_gen_mllg, mllg_gen);
  TH1D *hgen_mllg_EB_mu  = new TH1D("hgen_mllg_EB_mu", "hgen_mllg_EB_mu", bin_gen_mllg, mllg_gen);
  TH1D *hgen_mllg_EE_mu  = new TH1D("hgen_mllg_EE_mu", "hgen_mllg_EE_mu", bin_gen_mllg, mllg_gen);

  //--------------- pt llg ---------------------//

  const int bin_gen_ptllg = 8;
  float ptllg_gen[bin_gen_ptllg+1] = {0, 10, 20, 30, 40, 50, 70, 120, 1000};
  const int bin_rec_ptllg = 16;
  float ptllg_rec[bin_rec_ptllg+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 90, 120, 150, 1000};

  TH1D *hrec_ptllg_EB_ele = new TH1D("hrec_ptllg_EB_ele", "hrec_ptllg_EB_ele", bin_rec_ptllg, ptllg_rec);
  TH1D *hrec_ptllg_EE_ele = new TH1D("hrec_ptllg_EE_ele", "hrec_ptllg_EE_ele", bin_rec_ptllg, ptllg_rec);
  TH1D *hrec_ptllg_EB_mu = new TH1D("hrec_ptllg_EB_mu", "hrec_ptllg_EB_mu", bin_rec_ptllg, ptllg_rec);
  TH1D *hrec_ptllg_EE_mu = new TH1D("hrec_ptllg_EE_mu", "hrec_ptllg_EE_mu", bin_rec_ptllg, ptllg_rec);

  TH1D *hgen_ptllg_EB_ele = new TH1D("hgen_ptllg_EB_ele", "hgen_ptllg_EB_ele", bin_gen_ptllg, ptllg_gen);
  TH1D *hgen_ptllg_EE_ele = new TH1D("hgen_ptllg_EE_ele", "hgen_ptllg_EE_ele", bin_gen_ptllg, ptllg_gen);
  TH1D *hgen_ptllg_EB_mu  = new TH1D("hgen_ptllg_EB_mu", "hgen_ptllg_EB_mu", bin_gen_ptllg, ptllg_gen);
  TH1D *hgen_ptllg_EE_mu  = new TH1D("hgen_ptllg_EE_mu", "hgen_ptllg_EE_mu", bin_gen_ptllg, ptllg_gen);

  // --------------- jet multiplicity ------------//

  const int njet_recbin = 7;
  const int njet_genbin = 4;

  TH1D *hrec_ele_njet = new TH1D("hrec_ele_njet", "hrec_ele_njet", njet_recbin, 0, njet_recbin);
  TH1D *hgen_ele_njet = new TH1D("hgen_ele_njet", "hgen_ele_njet", njet_genbin, 0, njet_genbin);

  TH1D *hrec_mu_njet = new TH1D("hrec_mu_njet", "hrec_mu_njet", njet_recbin, 0, njet_recbin);
  TH1D *hgen_mu_njet = new TH1D("hgen_mu_njet", "hgen_mu_njet", njet_genbin, 0, njet_genbin);

  cout << "loop over mc tree" << endl;
  cout << "total events from tree: " << tmc->GetEntriesFast() << endl;

  int countEle = 0;
  int countMu = 0;
  int passEv = 0;

  float ev_jetrec[njet_recbin] = {0.};
  float ev_jetrec_mu[njet_recbin] = {0.};

  //Long64_t startEv = tmc->GetEntriesFast() * 0.8; //for amcatnlo sample
  Long64_t startEv = 0; //other samples
  if (sample.Contains("sherpa")) startEv = 0;
  //testing
  for (Long64_t iev = startEv; iev < tmc->GetEntriesFast(); iev++) {
    tmc->GetEntry(iev);
    if(tmc->GetEntry(iev)<=0) break;
    if (iev%20000 == 0) cout << "processing event: " << iev+1 << " th" << endl;

    float truth_Pt = genPhoEt;
    float truth_mllg = genXm;

    if (genPhoEt < 0. || genXm < 0.) continue;
    if ( fabs(genPhoEta)>2.5 ) continue;
    //if ( fabs(genlepEta[0]) > 2.4 || fabs(genlepEta[1]) > 2.4 ) continue;
    //if (genlepPt[0] < 25 || genlepPt[1] < 20) continue;
    if (genZm < 50.) continue;
    if (z_mass < 50. || z_charge != 0) continue;
    if (lept0_pt < 25. || lept1_pt < 20.) continue;
    if (fabs(lept0_eta) > 2.4 || fabs(lept1_eta) > 2.4) continue;
    if (gamma_pt < 0.) continue;
    //if (nselJet !=0 ) continue;

    passEv++;
    if (gamma_pt < 15.) gamma_SF = 1.;
    int ntruthJet = ngenjet; //for exclusive

    //electron 
    if ( trig_Ele23_Ele12 ==1 && leptType ==11) {

      float Xm = boss_mass;

      float wei = puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;
      float weigen = 1;
      if (gamma_pt >= 15.) weigen = puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;
      else weigen = puweigj_65nb*genWeight*lept0_RecoSF*lept1_RecoSF*lept0_SelSF*lept1_SelSF;

      if (isEB && gamma_ChIso<2.) {
        hrec_ele_EB->Fill(gamma_pt,wei);
        hgen_ele_EB->Fill(truth_Pt,weigen);

        if (gamma_pt >= 20.) {
          hrec_mllg_EB_ele->Fill(Xm,wei);
          hgen_mllg_EB_ele->Fill(truth_mllg,weigen);

          hrec_ptllg_EB_ele->Fill(boss_pt,wei);
          hgen_ptllg_EB_ele->Fill(genXpt,weigen);
        }
      }
      else if( isEE && gamma_ChIso<1.5)  {
        hrec_ele_EE->Fill(gamma_pt,wei);
        hgen_ele_EE->Fill(truth_Pt,weigen);

        if (gamma_pt >= 20.) {
	  hrec_mllg_EE_ele->Fill(Xm,wei);
	  hgen_mllg_EE_ele->Fill(truth_mllg,weigen);

	  hrec_ptllg_EE_ele->Fill(boss_pt,wei);
	  hgen_ptllg_EE_ele->Fill(genXpt,weigen);
        }
      }

      if ( gamma_pt >= 20. && ((isEB && gamma_ChIso<2.) || (isEE && gamma_ChIso<1.5)) ) {
	hgen_ele_njet->Fill(TMath::Min(ntruthJet,3),weigen);

        if (nselJet==0) ev_jetrec[0] += wei;
        if (nselJet==1) {
          if ( jetPt[0] < 50.) ev_jetrec[1] += wei;
          if (jetPt[0] > 50.) ev_jetrec[2] += wei;
        }
        if (nselJet==2) {
          if ( jetPt[1] < 50.) ev_jetrec[3] += wei;
          if (jetPt[1] > 50.) ev_jetrec[4] += wei;
        }
        if (nselJet > 2) {
          if ( jetPt[2] < 50.) ev_jetrec[5] += wei;
          if (jetPt[2] > 50.) ev_jetrec[6] += wei;
        }
      }
    }

    //muon
    if (trig_Mu17_Mu8==1 && leptType==13 && pair_dPhi>70.) {
      float Xm = boss_mass;

      float wei = puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;
      float weigen = 1;
      if (gamma_pt >= 15.) weigen = puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF*gamma_SF*gamma_CSEV_SF;
      else weigen = puweigj_65nb*genWeight*lept0_SelSF*lept1_SelSF;
      if (isEB && gamma_ChIso<2.) {
        hrec_mu_EB->Fill(gamma_pt,wei);
        hgen_mu_EB->Fill(truth_Pt,weigen);

        if (gamma_pt >= 20.) {
	  hrec_mllg_EB_mu->Fill(Xm,wei);
	  hgen_mllg_EB_mu->Fill(truth_mllg,weigen);

	  hrec_ptllg_EB_mu->Fill(boss_pt,wei);
	  hgen_ptllg_EB_mu->Fill(genXpt,weigen);
        }
      }
      else if( isEE && gamma_ChIso<1.5)  {
        hrec_mu_EE->Fill(gamma_pt,wei);
        hgen_mu_EE->Fill(truth_Pt,weigen);

        if (gamma_pt >= 20.) {
	  hrec_mllg_EE_mu->Fill(Xm,wei);
	  hgen_mllg_EE_mu->Fill(truth_mllg,weigen);

	  hrec_ptllg_EE_mu->Fill(boss_pt,wei);
	  hgen_ptllg_EE_mu->Fill(genXpt,weigen);
        }
      }

      if ( gamma_pt >= 20. && ((isEB && gamma_ChIso<2.) || (isEE && gamma_ChIso<1.5)) ) {
        hgen_mu_njet->Fill(TMath::Min(ntruthJet,3),weigen);

        if (nselJet==0) ev_jetrec_mu[0] += wei;
	if (nselJet==1) {
          if ( jetPt[0] < 50.) ev_jetrec_mu[1] += wei;
          if (jetPt[0] > 50.) ev_jetrec_mu[2] += wei;
        }
        if (nselJet==2) {
          if ( jetPt[1] < 50.) ev_jetrec_mu[3] += wei;
          if (jetPt[1] > 50.) ev_jetrec_mu[4] += wei;
	}
        if (nselJet > 2) {
          if ( jetPt[2] < 50.) ev_jetrec_mu[5] += wei;
          if (jetPt[2] > 50.) ev_jetrec_mu[6] += wei;
        }
      }

    }
  }

  //jet multiplicity at reco level
  for (int i = 0; i < njet_recbin; i++) {
    hrec_ele_njet->SetBinContent(i+1, ev_jetrec[i]);
    hrec_ele_njet->SetBinError(i+1, TMath::Sqrt(ev_jetrec[i]));

    hrec_mu_njet->SetBinContent(i+1, ev_jetrec_mu[i]);
    hrec_mu_njet->SetBinError(i+1, TMath::Sqrt(ev_jetrec_mu[i]));

    cout << "rec jet of ele: " << ev_jetrec[i] << endl;
  }


  cout << "total ele channel: " << countEle << endl;
  cout << "total muon channel: " << countMu << endl;
  cout << "total event passing: " << passEv << endl;

  //input from data
  TH1D *hda_ele_EB = new TH1D("hda_ele_EB", "hda_ele_EB", bin_rec, pt_rec);
  TH1D *hda_ele_EE = new TH1D("hda_ele_EE", "hda_ele_EE", bin_rec, pt_rec);

  TH1D *hda_mu_EB = new TH1D("hda_mu_EB", "hda_mu_EB", bin_rec, pt_rec);
  TH1D *hda_mu_EE = new TH1D("hda_mu_EE", "hda_mu_EE", bin_rec, pt_rec);

  float yield_da_ele_EB[bin_rec] = {0.};
  float yield_da_ele_EE[bin_rec] = {0.};
  float yield_da_mu_EB[bin_rec] = {0.};
  float yield_da_mu_EE[bin_rec] = {0.};
  float erryield_da_ele_EB[bin_rec] = {0.};
  float erryield_da_ele_EE[bin_rec] = {0.};
  float erryield_da_mu_EB[bin_rec] = {0.};
  float erryield_da_mu_EE[bin_rec] = {0.};

  TString sufix = "_SB_EB7to13_EE6to14_biasUp";
  ifstream Yieldfile_ptg_ele("fittingYield/yield_PhoPt_ele" + sufix + ".txt");
  ifstream Yieldfile_ptg_mu("fittingYield/yield_PhoPt_mu" + sufix + ".txt");

  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_ele >> yield_da_ele_EB[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_ele >> erryield_da_ele_EB[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_ele >> yield_da_ele_EE[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_ele >> erryield_da_ele_EE[ibin];

  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_mu >> yield_da_mu_EB[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_mu >> erryield_da_mu_EB[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_mu >> yield_da_mu_EE[ibin];
  for (int ibin=0; ibin < bin_rec; ibin++) Yieldfile_ptg_mu >> erryield_da_mu_EE[ibin];


  //filling histogram for data
  for (int i = 1; i <= bin_rec; i++) {
    hda_ele_EB->SetBinContent(i, yield_da_ele_EB[i-1]);
    hda_ele_EE->SetBinContent(i, yield_da_ele_EE[i-1]);

    hda_ele_EB->SetBinError(i, erryield_da_ele_EB[i-1]);
    hda_ele_EE->SetBinError(i, erryield_da_ele_EE[i-1]);

    hda_mu_EB->SetBinContent(i, yield_da_mu_EB[i-1]);
    hda_mu_EE->SetBinContent(i, yield_da_mu_EE[i-1]);

    hda_mu_EB->SetBinError(i, erryield_da_mu_EB[i-1]);
    hda_mu_EE->SetBinError(i, erryield_da_mu_EE[i-1]);

  }

  float yield_da_mllg_EB_ele[bin_rec_mllg] = {0.};
  float yield_da_mllg_EE_ele[bin_rec_mllg] = {0.};
  float yield_da_mllg_EB_mu[bin_rec_mllg] = {0.};
  float yield_da_mllg_EE_mu[bin_rec_mllg] = {0.};

  float erryield_da_mllg_EB_ele[bin_rec_mllg] = {0.};
  float erryield_da_mllg_EE_ele[bin_rec_mllg] = {0.};
  float erryield_da_mllg_EB_mu[bin_rec_mllg] = {0.};
  float erryield_da_mllg_EE_mu[bin_rec_mllg] = {0.};

  TH1D *hda_mllg_EB_ele = new TH1D("hda_mllg_EB_ele", "hda_mllg_EB_ele", bin_rec_mllg, mllg_rec);
  TH1D *hda_mllg_EE_ele = new TH1D("hda_mllg_EE_ele", "hda_mllg_EE_ele", bin_rec_mllg, mllg_rec);
  TH1D *hda_mllg_EB_mu = new TH1D("hda_mllg_EB_mu", "hda_mllg_EB_mu", bin_rec_mllg, mllg_rec);
  TH1D *hda_mllg_EE_mu = new TH1D("hda_mllg_EE_mu", "hda_mllg_EE_mu", bin_rec_mllg, mllg_rec);

  ifstream Yieldfile_mllg_ele("fittingYield/yield_Mllg_ele" + sufix + ".txt");
  ifstream Yieldfile_mllg_mu("fittingYield/yield_Mllg_mu" + sufix + ".txt");

  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_ele >> yield_da_mllg_EB_ele[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_ele >> erryield_da_mllg_EB_ele[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_ele >> yield_da_mllg_EE_ele[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_ele >> erryield_da_mllg_EE_ele[ibin];

  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_mu >> yield_da_mllg_EB_mu[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_mu >> erryield_da_mllg_EB_mu[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_mu >> yield_da_mllg_EE_mu[ibin];
  for (int ibin=0; ibin < bin_rec_mllg; ibin++) Yieldfile_mllg_mu >> erryield_da_mllg_EE_mu[ibin];


  for (int i = 1; i <= bin_rec_mllg; i++) {
    hda_mllg_EB_ele->SetBinContent(i, yield_da_mllg_EB_ele[i-1]);
    hda_mllg_EE_ele->SetBinContent(i, yield_da_mllg_EE_ele[i-1]);
    hda_mllg_EB_mu->SetBinContent(i, yield_da_mllg_EB_mu[i-1]);
    hda_mllg_EE_mu->SetBinContent(i, yield_da_mllg_EE_mu[i-1]);

    hda_mllg_EB_ele->SetBinError(i, erryield_da_mllg_EB_ele[i-1]);
    hda_mllg_EE_ele->SetBinError(i, erryield_da_mllg_EE_ele[i-1]);
    hda_mllg_EB_mu->SetBinError(i, erryield_da_mllg_EB_mu[i-1]);
    hda_mllg_EE_mu->SetBinError(i, erryield_da_mllg_EE_mu[i-1]);
  }

  
  //for Ptllg
  float yield_da_ptllg_EB_ele[bin_rec_ptllg] = {0.};
  float yield_da_ptllg_EE_ele[bin_rec_ptllg] = {0.};
  float yield_da_ptllg_EB_mu[bin_rec_ptllg] = {0.};
  float yield_da_ptllg_EE_mu[bin_rec_ptllg] = {0.};

  float erryield_da_ptllg_EB_ele[bin_rec_ptllg] = {0.};
  float erryield_da_ptllg_EE_ele[bin_rec_ptllg] = {0.};
  float erryield_da_ptllg_EB_mu[bin_rec_ptllg] = {0.};
  float erryield_da_ptllg_EE_mu[bin_rec_ptllg] = {0.};

  TH1D *hda_ptllg_EB_ele = new TH1D("hda_ptllg_EB_ele", "hda_ptllg_EB_ele", bin_rec_ptllg, ptllg_rec);
  TH1D *hda_ptllg_EE_ele = new TH1D("hda_ptllg_EE_ele", "hda_ptllg_EE_ele", bin_rec_ptllg, ptllg_rec);
  TH1D *hda_ptllg_EB_mu = new TH1D("hda_ptllg_EB_mu", "hda_ptllg_EB_mu", bin_rec_ptllg, ptllg_rec);
  TH1D *hda_ptllg_EE_mu = new TH1D("hda_ptllg_EE_mu", "hda_ptllg_EE_mu", bin_rec_ptllg, ptllg_rec);

  ifstream Yieldfile_ptllg_ele("fittingYield/yield_Ptllg_ele" + sufix + ".txt");
  ifstream Yieldfile_ptllg_mu("fittingYield/yield_Ptllg_mu" + sufix + ".txt");

  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_ele >> yield_da_ptllg_EB_ele[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_ele >> erryield_da_ptllg_EB_ele[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_ele >> yield_da_ptllg_EE_ele[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_ele >> erryield_da_ptllg_EE_ele[ibin];

  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_mu >> yield_da_ptllg_EB_mu[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_mu >> erryield_da_ptllg_EB_mu[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_mu >> yield_da_ptllg_EE_mu[ibin];
  for (int ibin=0; ibin < bin_rec_ptllg; ibin++) Yieldfile_ptllg_mu >> erryield_da_ptllg_EE_mu[ibin];


  for (int i = 1; i <= bin_rec_ptllg; i++) {
    hda_ptllg_EB_ele->SetBinContent(i, yield_da_ptllg_EB_ele[i-1]);
    hda_ptllg_EE_ele->SetBinContent(i, yield_da_ptllg_EE_ele[i-1]);
    hda_ptllg_EB_mu->SetBinContent(i, yield_da_ptllg_EB_mu[i-1]);
    hda_ptllg_EE_mu->SetBinContent(i, yield_da_ptllg_EE_mu[i-1]);

    hda_ptllg_EB_ele->SetBinError(i, erryield_da_ptllg_EB_ele[i-1]);
    hda_ptllg_EE_ele->SetBinError(i, erryield_da_ptllg_EE_ele[i-1]);
    hda_ptllg_EB_mu->SetBinError(i, erryield_da_ptllg_EB_mu[i-1]);
    hda_ptllg_EE_mu->SetBinError(i, erryield_da_ptllg_EE_mu[i-1]);
  }
  
  
  //for exclusive , no need multiplicity

  //-------------------------------------------------------//
  //yield for jet multiplicity
  TH1D *hda_njet_ele = new TH1D("hda_njet_ele", "hda_njet_ele", njet_recbin, 0, njet_recbin);
  TH1D *hda_njet_mu = new TH1D("hda_njet_mu", "hda_njet_mu", njet_recbin, 0, njet_recbin);

  float yield_da_njet_ele[njet_recbin] = {0.};
  float yield_da_njet_mu[njet_recbin] = {0.};
  float erryield_da_njet_ele[njet_recbin] = {0.};
  float erryield_da_njet_mu[njet_recbin] = {0.};

  ifstream Yieldfile_excl_njet_ele("fittingYield/yield_excl_njet_ele" + sufix + ".txt");
  ifstream Yieldfile_excl_njet_mu("fittingYield/yield_excl_njet_mu" + sufix + ".txt");

  for (int ibin=0; ibin < njet_recbin; ibin++) Yieldfile_excl_njet_ele >> yield_da_njet_ele[ibin];
  for (int ibin=0; ibin < njet_recbin; ibin++) Yieldfile_excl_njet_ele >> erryield_da_njet_ele[ibin];

  for (int ibin=0; ibin < njet_recbin; ibin++) Yieldfile_excl_njet_mu >> yield_da_njet_mu[ibin];
  for (int ibin=0; ibin < njet_recbin; ibin++) Yieldfile_excl_njet_mu >> erryield_da_njet_mu[ibin];

  for (int i = 1; i<= njet_recbin; i++){
    hda_njet_ele->SetBinContent(i, yield_da_njet_ele[i-1]);
    hda_njet_mu->SetBinContent(i, yield_da_njet_mu[i-1]);

    hda_njet_ele->SetBinError(i, erryield_da_njet_ele[i-1]);
    hda_njet_mu->SetBinError(i, erryield_da_njet_mu[i-1]);
  }

  
  //for output 
  TString outname = "inputUnfold/SB_EE7to13_EE6to14/main/inputUnfold_isPVGood";
  if (sample.Contains("amcnlo")) outname += "_Madgraph_data";
  else outname += "_sherpa";

  if (sufix.Contains("biasUp")) outname += "_biasUp";
  if (sufix.Contains("biasDn")) outname += "_biasDn";
  //outname += "_Mllg70_StatNJet.root";
  //outname += "_ttg_vv_bkg.root";
  outname += ".root";

  //other background 
  TFile *fbkg_EBpho_ele = new TFile("Bkg/TTG_VV_bkg_EBpho_eleChan.root", "read");
  TFile *fbkg_EEpho_ele = new TFile("Bkg/TTG_VV_bkg_EEpho_eleChan.root");
  TFile *fbkg_EBpho_mu = new TFile("Bkg/TTG_VV_bkg_EBpho_muChan.root");
  TFile *fbkg_EEpho_mu = new TFile("Bkg/TTG_VV_bkg_EEpho_muChan.root");

  TH1F *hbkg_ttg_vv_ptg_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_mllg_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_ptllg_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_ttg_vv_ptllg");

  TH1F *hbkg_ttg_vv_ptg_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_mllg_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_ptllg_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_ttg_vv_ptllg");

  TH1F *hbkg_ttg_vv_ptg_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_mllg_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_ptllg_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_ttg_vv_ptllg");

  TH1F *hbkg_ttg_vv_ptg_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_ttg_vv_ptg");
  TH1F *hbkg_ttg_vv_mllg_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_ttg_vv_mllg");
  TH1F *hbkg_ttg_vv_ptllg_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_ttg_vv_ptllg");


  TFile *fbkg_njet = new TFile("Bkg/TTG_VV_bkg_njet.root", "read");
  TH1F *hbkg_ttg_vv_njet_ele = (TH1F*)fbkg_njet->Get("hbkg_ttg_vv_njet_ele");
  TH1F *hbkg_ttg_vv_njet_mu = (TH1F*)fbkg_njet->Get("hbkg_ttg_vv_njet_mu");

  /* // rum in cms01, before July 15
  //inclusive
  TFile *fbkg_EBpho_ele = new TFile("~/analysis/Zg_SM/2016data/Diboson_ttb_bkg/diboson_bgk_EBpho_eleChan.root", "read");
  TFile *fbkg_EEpho_ele = new TFile("~/analysis/Zg_SM/2016data/Diboson_ttb_bkg/diboson_bgk_EEpho_eleChan.root", "read");
  TFile *fbkg_EBpho_mu = new TFile("~/analysis/Zg_SM/2016data/Diboson_ttb_bkg/diboson_bgk_EBpho_muChan.root", "read");
  TFile *fbkg_EEpho_mu = new TFile("~/analysis/Zg_SM/2016data/Diboson_ttb_bkg/diboson_bgk_EEpho_muChan.root", "read");


  TH1F *hbkg_diboson_pho_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_diboson_pho");
  TH1F *hbkg_diboson_mllg_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_diboson_mllg");
  TH1F *hbkg_ttb_pho_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_ttb_pho");
  TH1F *hbkg_ttb_mllg_EB_ele = (TH1F*) fbkg_EBpho_ele->Get("hbkg_ttb_mllg");

  TH1F *hbkg_diboson_pho_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_diboson_pho");
  TH1F *hbkg_diboson_mllg_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_diboson_mllg");
  TH1F *hbkg_ttb_pho_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_ttb_pho");
  TH1F *hbkg_ttb_mllg_EE_ele = (TH1F*) fbkg_EEpho_ele->Get("hbkg_ttb_mllg");

  TH1F *hbkg_diboson_pho_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_diboson_pho");
  TH1F *hbkg_diboson_mllg_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_diboson_mllg");
  TH1F *hbkg_ttb_pho_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_ttb_pho");
  TH1F *hbkg_ttb_mllg_EB_mu = (TH1F*) fbkg_EBpho_mu->Get("hbkg_ttb_mllg");

  TH1F *hbkg_diboson_pho_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_diboson_pho");
  TH1F *hbkg_diboson_mllg_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_diboson_mllg");
  TH1F *hbkg_ttb_pho_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_ttb_pho");
  TH1F *hbkg_ttb_mllg_EE_mu = (TH1F*) fbkg_EEpho_mu->Get("hbkg_ttb_mllg");

  TFile *fbkg_njet = new TFile("~/analysis/Zg_SM/2016data/yield_njet/diboson_ttb_bkg_njet.root", "read");
  TH1F *hbkg_diboson_njet_ele = (TH1F*) fbkg_njet->Get("hbkg_diboson_njet_ele");
  TH1F *hbkg_ttb_njet_ele = (TH1F*) fbkg_njet->Get("hbkg_ttb_njet_ele");
  TH1F *hbkg_diboson_njet_mu = (TH1F*) fbkg_njet->Get("hbkg_diboson_njet_mu");
  TH1F *hbkg_ttb_njet_mu = (TH1F*) fbkg_njet->Get("hbkg_ttb_njet_mu");
  */

  
  TFile *outf = new TFile(outname, "recreate");
  outf->cd();
  
  //reco MC for testing
  hrec_ele_EB->Write();
  hrec_ele_EE->Write();
  hrec_mu_EB->Write();
  hrec_mu_EE->Write();
  hrec_mllg_EB_ele->Write();
  hrec_mllg_EE_ele->Write();
  hrec_mllg_EB_mu->Write();
  hrec_mllg_EE_mu->Write();

  hrec_ptllg_EB_ele->Write();
  hrec_ptllg_EE_ele->Write();
  hrec_ptllg_EB_mu->Write();
  hrec_ptllg_EE_mu->Write();
  hrec_ele_njet->SetName("hrec_njet_ele");
  hrec_mu_njet->SetName("hrec_njet_mu");
  hrec_ele_njet->Write();
  hrec_mu_njet->Write();

  //gen mc for testing
  hgen_ele_EB->Write();
  hgen_ele_EE->Write();
  hgen_mu_EB->Write();
  hgen_mu_EE->Write();
  hgen_mllg_EB_ele->Write();
  hgen_mllg_EE_ele->Write();
  hgen_mllg_EB_mu->Write();
  hgen_mllg_EE_mu->Write();

  hgen_ptllg_EB_ele->Write();
  hgen_ptllg_EE_ele->Write();
  hgen_ptllg_EB_mu->Write();
  hgen_ptllg_EE_mu->Write();

  hgen_ele_njet->SetName("hgen_njet_ele");
  hgen_mu_njet->SetName("hgen_njet_mu");
  hgen_ele_njet->Write();
  hgen_mu_njet->Write();

  
  //data
  hda_ele_EB->Write();
  hda_ele_EE->Write();
  hda_mu_EB->Write();
  hda_mu_EE->Write();

  hda_mllg_EB_ele->Write();
  hda_mllg_EE_ele->Write();
  hda_mllg_EB_mu->Write();
  hda_mllg_EE_mu->Write();

  hda_ptllg_EB_ele->Write();
  hda_ptllg_EE_ele->Write();
  hda_ptllg_EB_mu->Write();
  hda_ptllg_EE_mu->Write();

  hda_njet_ele->Write();
  hda_njet_mu->Write();

  hbkg_ttg_vv_ptg_EB_ele->Write("hbkg_ttg_vv_ptg_EB_ele");
  hbkg_ttg_vv_mllg_EB_ele->Write("hbkg_ttg_vv_mllg_EB_ele");
  hbkg_ttg_vv_ptllg_EB_ele->Write("hbkg_ttg_vv_ptllg_EB_ele");

  hbkg_ttg_vv_ptg_EE_ele->Write("hbkg_ttg_vv_ptg_EE_ele");
  hbkg_ttg_vv_mllg_EE_ele->Write("hbkg_ttg_vv_mllg_EE_ele");
  hbkg_ttg_vv_ptllg_EE_ele->Write("hbkg_ttg_vv_ptllg_EE_ele");

  hbkg_ttg_vv_ptg_EB_mu->Write("hbkg_ttg_vv_ptg_EB_mu");
  hbkg_ttg_vv_mllg_EB_mu->Write("hbkg_ttg_vv_mllg_EB_mu");
  hbkg_ttg_vv_ptllg_EB_mu->Write("hbkg_ttg_vv_ptllg_EB_mu");

  hbkg_ttg_vv_ptg_EE_mu->Write("hbkg_ttg_vv_ptg_EE_mu");
  hbkg_ttg_vv_mllg_EE_mu->Write("hbkg_ttg_vv_mllg_EE_mu");
  hbkg_ttg_vv_ptllg_EE_mu->Write("hbkg_ttg_vv_ptllg_EE_mu");

  hbkg_ttg_vv_njet_ele->Write("hbkg_ttg_vv_njet_ele");
  hbkg_ttg_vv_njet_mu->Write("hbkg_ttg_vv_njet_mu");
  

  //background
  /* //before July 20
  hbkg_diboson_pho_EB_ele->Write("hbkg_diboson_pho_EB_ele");
  hbkg_diboson_pho_EE_ele->Write("hbkg_diboson_pho_EE_ele");
  hbkg_diboson_mllg_EB_ele->Write("hbkg_diboson_mllg_EB_ele");
  hbkg_diboson_mllg_EE_ele->Write("hbkg_diboson_mllg_EE_ele");

  hbkg_ttb_pho_EB_ele->Write("hbkg_ttb_pho_EB_ele");
  hbkg_ttb_pho_EE_ele->Write("hbkg_ttb_pho_EE_ele");
  hbkg_ttb_mllg_EB_ele->Write("hbkg_ttb_mllg_EB_ele");
  hbkg_ttb_mllg_EE_ele->Write("hbkg_ttb_mllg_EE_ele");

  hbkg_diboson_pho_EB_mu->Write("hbkg_diboson_pho_EB_mu");
  hbkg_diboson_pho_EE_mu->Write("hbkg_diboson_pho_EE_mu");
  hbkg_diboson_mllg_EB_mu->Write("hbkg_diboson_mllg_EB_mu");
  hbkg_diboson_mllg_EE_mu->Write("hbkg_diboson_mllg_EE_mu");

  hbkg_ttb_pho_EB_mu->Write("hbkg_ttb_pho_EB_mu");
  hbkg_ttb_pho_EE_mu->Write("hbkg_ttb_pho_EE_mu");
  hbkg_ttb_mllg_EB_mu->Write("hbkg_ttb_mllg_EB_mu");
  hbkg_ttb_mllg_EE_mu->Write("hbkg_ttb_mllg_EE_mu");


  hbkg_diboson_njet_ele->Write("hbkg_diboson_njet_ele");
  hbkg_ttb_njet_ele->Write("hbkg_ttb_njet_ele");

  hbkg_diboson_njet_mu->Write("hbkg_diboson_njet_mu");
  hbkg_ttb_njet_mu->Write("hbkg_ttb_njet_mu");
  */

  outf->Write();
  outf->Close();
  fmc->Close();

  cout << "done" << endl;

}
