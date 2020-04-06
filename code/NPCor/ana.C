#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLorentzVector.h"

#include <iostream>
using std::cout;
#endif

using namespace std;
using namespace reco;

float deltaPhi(float phi1, float phi2) {

  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return dPhi;

}

float deltaR(float eta1, float phi1, float eta2, float phi2) {

  float dEta = eta1 - eta2;
  float dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2*TMath::Pi();
  if (dPhi <= -TMath::Pi()) dPhi += 2* TMath::Pi();

  return sqrt(pow(dEta,2) + pow (dPhi,2));

}

float getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
                   reco::GenParticleCollection::const_iterator thisPart,
                   float dRMax, bool removeMu, bool removeNu) {

  // Returns Et sum
  float etSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
      continue;

    int pdgCode = abs(p->pdgId());

    // skip muons/neutrinos, if requested
    if (removeMu && pdgCode == 13) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;

    // must be within deltaR cone
    //float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    float dR = deltaR(thisPart->eta(), thisPart->phi(), p->eta(), p->phi());
    if (dR > dRMax) continue;

    etSum += p->et();
  }

  return etSum;
}


void ana(TString filename = "/data3/ggNtuples/Zg_summer16/Zg_Madgraph_0123J_5f_MoreStat.root") {

  TString outname = "Zg_Madgraph_0123J_5f_7M_QCD.root";

  TFile *fout = new TFile(outname, "recreate");

  TTree *outtree = new TTree("outtree", "output tree");
  float pho_pt;
  float pho_eta;
  float pho_phi;
  float ele1_pt;
  float ele2_pt;
  float ele1_eta;
  float ele2_eta;
  float ele1_phi;
  float ele2_phi;
  float mu1_eta;
  float mu2_eta;
  float mu1_phi;
  float mu2_phi;
  float mu1_pt;
  float mu2_pt;
  float genWeight;
  float boss_mass;
  float boss_pt;
  float z_mass;
  float z_pt;
  float njet;

  outtree->Branch("pho_pt",    &pho_pt);
  outtree->Branch("pho_eta",   &pho_eta);
  outtree->Branch("pho_phi",   &pho_phi);
  outtree->Branch("ele1_pt",    &ele1_pt);
  outtree->Branch("ele1_eta",   &ele1_eta);
  outtree->Branch("ele1_phi",   &ele1_phi);
  outtree->Branch("ele2_pt",    &ele2_pt);
  outtree->Branch("ele2_eta",   &ele2_eta);
  outtree->Branch("ele2_phi",   &ele2_phi);
  outtree->Branch("mu1_pt",    &mu1_pt);
  outtree->Branch("mu1_eta",   &mu1_eta);
  outtree->Branch("mu1_phi",   &mu1_phi);
  outtree->Branch("mu2_pt",    &mu2_pt);
  outtree->Branch("mu2_eta",   &mu2_eta);
  outtree->Branch("mu2_phi",   &mu2_phi);
  outtree->Branch("boss_mass", &boss_mass);
  outtree->Branch("boss_pt",   &boss_pt);
  outtree->Branch("z_mass",    &z_mass);
  outtree->Branch("z_pt",      &z_pt);
  outtree->Branch("njet",      &njet);
  outtree->Branch("genWeight", &genWeight);

  int nbin_pt = 10;
  float ptbin[] = {0, 10, 20, 25, 30, 35, 45, 60, 80, 120, 1000};
  TH1F *hPhoPt = new TH1F("hPhoPt", "hPhoPt", nbin_pt, ptbin);

  int nbin_m = 10;
  float mllg[] = {70, 88, 95, 110, 135, 170, 210, 270, 350, 470, 1000};
  TH1F *hMllg = new TH1F("hMllg", "hMllg", nbin_m, mllg);
  TH1D *htotwei = new TH1D("htotwei", "htotwei", 2, 0, 2);

  TFile file(filename);

  fwlite::Event ev(&file);

  double totWei = 0.;
  int ipho = 0;

  Long64_t nEvent = 0;
  for( ev.toBegin(); ! ev.atEnd(); ++ev) {
    //if (nEvent > 20000) continue;
    if ((nEvent%10000) == 0) cout << "processing event: " << nEvent << " th" << endl;
    nEvent++;

    edm::EventBase const & event = ev;
    edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
    event.getByLabel(std::string("genParticles"), genParticlesHandle);

    TLorentzVector pho, lep1, lep2;
    int nele = 0, nmu = 0, npho = 0;
    pho_pt = 0.;
    ele1_pt = 0., ele2_pt = 0.;
    mu1_pt = 0., mu2_pt = 0.;
    genWeight = 0.;
    ele1_eta = -99., ele1_phi = -99.;
    ele2_eta = -99., ele2_phi = -99.;
    mu1_eta = -99., mu1_phi = -99.;
    mu2_eta = -99., mu2_phi = -99.;
    pho_eta = -99., pho_phi = -99.;
    njet = 0;
    float genPhoIso = 999.;
    vector<float> vec_pho_pt, vec_pho_eta, vec_pho_phi;
    vector<float> vec_ele_pt, vec_ele_eta, vec_ele_phi;
    vector<float> vec_mu_pt, vec_mu_eta, vec_mu_phi;
    vec_pho_pt.clear();
    vec_pho_eta.clear();
    vec_pho_phi.clear();

    vec_ele_pt.clear();
    vec_ele_eta.clear();
    vec_ele_phi.clear();

    vec_mu_pt.clear();
    vec_mu_eta.clear();
    vec_mu_phi.clear();


    for(vector<reco::GenParticle>::const_iterator igen=genParticlesHandle->begin(); igen!=genParticlesHandle->end(); ++igen){
      int status = igen->status();
      int pdgId = igen->pdgId();
      int isPrompt = igen->isPromptFinalState();
      int isHadProcess = igen->fromHardProcessFinalState();
      //cout << "status: " << status << endl;
      if (status != 1) continue;

      //cout << "loop on gen particle level" << endl;
      if (abs(pdgId) == 22 && igen->pt()>10) {
        genPhoIso = getGenCalIso(genParticlesHandle, igen, 0.4, false, false);
        if (genPhoIso > 5.) continue;
        if (fabs(igen->eta()) > 2.5) continue;
        if (igen->pt() < 15.) continue;
        pho.SetPtEtaPhiM(igen->pt(), igen->eta(), igen->phi(), 0.0);
	vec_pho_pt.push_back(igen->pt());
        vec_pho_eta.push_back(igen->eta());
        vec_pho_phi.push_back(igen->phi());

        ipho++;
      }

      if (abs(pdgId) == 11 && igen->pt()>20) {
        if ( fabs(igen->eta()) > 2.5) continue;
	vec_ele_pt.push_back(igen->pt());
        vec_ele_eta.push_back(igen->eta());
        vec_ele_phi.push_back(igen->phi());

      }
      if (abs(pdgId) == 13 && igen->pt()>20) {
	if ( fabs(igen->eta()) > 2.5) continue;
	vec_mu_pt.push_back(igen->pt());
        vec_mu_eta.push_back(igen->eta());
        vec_mu_phi.push_back(igen->phi());
      }
    }

    //sort in pt
    if ( vec_pho_pt.size() < 1) continue;

    //cout << "npho before sorting: "<< vec_pho_pt.size() << "\t nele: " << vec_ele_pt.size() << endl;
    float maxpt = 0;
    for (unsigned int i = 0; i < vec_pho_pt.size(); i++) {
      //cout << "pho pt: " << vec_pho_pt[i] << endl;
      if (maxpt < vec_pho_pt[i]) {
        maxpt = vec_pho_pt[i];
        pho_pt = vec_pho_pt[i];
        pho.SetPtEtaPhiM(vec_pho_pt[i], vec_pho_eta[i], vec_pho_phi[i], 0.);
      }
    }
    
    if ( vec_ele_pt.size() < 2 && vec_mu_pt.size() < 2) continue;

    float max_ele1_pt = 0;
    float max_ele2_pt = 0;
    unsigned int index_ele1 = -1;
    for (unsigned int i = 0; i < vec_ele_pt.size(); i++) {
      //cout << "ele pt: " << vec_ele_pt[i] << endl;
      if (max_ele1_pt < vec_ele_pt[i]) {
	max_ele1_pt = vec_ele_pt[i];
	lep1.SetPtEtaPhiM(vec_ele_pt[i], vec_ele_eta[i], vec_ele_phi[i], 0.511*0.001);
        index_ele1 = i;
      }
    }

    //cout << "sort 2nd ele" << endl;
    for (unsigned int i = 0; i < vec_ele_pt.size(); i++) {
      if (i == index_ele1) continue;
      if ( max_ele2_pt < vec_ele_pt[i]) {
        max_ele2_pt = vec_ele_pt[i];
	lep2.SetPtEtaPhiM(vec_ele_pt[i], vec_ele_eta[i], vec_ele_phi[i], 0.511*0.001);
      }
    }

    //cout << "ele1: " << max_ele1_pt << "\t ele2: " << max_ele2_pt << "\t pho: " << pho_pt << endl;

    float max_mu1_pt = 0;
    float max_mu2_pt = 0;
    unsigned int index_mu1 = -1;
    for (unsigned int i = 0; i < vec_mu_pt.size(); i++) {
      if (max_mu1_pt < vec_mu_pt[i]) {
        max_mu1_pt = vec_mu_pt[i];
	lep1.SetPtEtaPhiM(vec_mu_pt[i], vec_mu_eta[i], vec_mu_phi[i], 104*0.001);
        index_mu1 = i;
      }
    }

    for (unsigned int i = 0; i < vec_mu_pt.size(); i++) {
      if (i == index_mu1) continue;
      if ( max_mu2_pt < vec_mu_pt[i]) {
        max_mu2_pt = vec_mu_pt[i];
	lep2.SetPtEtaPhiM(vec_mu_pt[i], vec_mu_eta[i], vec_mu_phi[i], 104*0.001);
      }
    }


    if (deltaR(pho.Eta(), pho.Phi(), lep1.Eta(), lep1.Phi()) < 0.7 ) continue;
    if (deltaR(pho.Eta(), pho.Phi(), lep2.Eta(), lep2.Phi()) < 0.7 ) continue;

    edm::Handle<GenEventInfoProduct> genEventInfoHandle;
    event.getByLabel(std::string("generator"), genEventInfoHandle);

    float genWeight_ = 1;
    if (genEventInfoHandle.isValid())
      genWeight_ = genEventInfoHandle->weight();

    if (genWeight_ > 0) genWeight = 1;
    else genWeight = -1;
    totWei += genWeight;

    //get genjet
    edm::Handle<vector<reco::GenJet> > genJetHandle;
    event.getByLabel(std::string("ak4GenJets"), genJetHandle);

    int ngenjet = 0;
    for(vector<reco::GenJet>::const_iterator igen=genJetHandle->begin(); igen!=genJetHandle->end(); ++igen){
      int status = igen->status();
      int pdgId = abs(igen->pdgId());
      if (status != 1) continue;
      if ( (pdgId >=1 && pdgId <=6) || pdgId == 21) {//quarks n gluons
	if ( igen->pt() < 30.) continue;
        if (fabs(igen->eta()) > 2.4) continue;
        if (deltaR(pho.Eta(), pho.Phi(), igen->eta(), igen->phi()) < 0.4) continue;
	if (deltaR(lep1.Eta(), lep1.Phi(), igen->eta(), igen->phi()) < 0.4) continue;
        if (deltaR(lep2.Eta(), lep2.Phi(), igen->eta(), igen->phi()) < 0.4) continue;

        ngenjet++;
      }
    }
    //cout << "ngenjet: " << ngenjet << endl;
    //if (ngenjet > 0) continue;
    njet = ngenjet;

    //cout << "ele1 max: " << max_ele1_pt << "\t from lorentz: " << lep1.Pt() << endl;
    //cout << "pho from lorentz: " << pho.Pt() << "\t from other: " << pho_pt << endl;
    if (pho.Pt()>20. && ((max_ele1_pt>25. && max_ele2_pt>20.) || (max_mu1_pt>25. && max_mu2_pt>20.)) ) {

      TLorentzVector X = pho + lep1 + lep2;
      TLorentzVector Z = lep1 + lep2;
      if (Z.M() < 50.) continue;
      if (fabs(pho.Eta()) > 2.5 || fabs(lep1.Eta())>2.4 || fabs(lep2.Eta())>2.4) continue;

      boss_mass = X.M();
      boss_pt = X.Pt();
      z_mass = (lep1+lep2).M();
      z_pt = (lep1+lep2).Pt();
      hPhoPt->Fill(pho.Pt(), genWeight);
      hMllg->Fill(X.M(), genWeight);

      pho_pt = pho.Pt();
      pho_eta = pho.Eta();
      pho_phi = pho.Phi();
      ele1_pt = lep1.Pt();
      ele1_eta = lep1.Eta();
      ele1_phi = lep1.Phi();
      ele2_pt = lep2.Pt();
      ele2_eta = lep2.Eta();
      ele2_phi = lep2.Phi();


      outtree->Fill();
    }
  }
  cout << "done event loop" << endl;

  htotwei->SetBinContent(1,totWei);

  //float xs = 47.33e3; //Zg_nlo_noqed
  //float xs = 47.75e3; //Zg_nlo_nofsr
  float xs = 37.92e3; //Zg_lo_noqcd
  //float xs = 47.42e3; // noqed_withFsr

  float scale = totWei/xs;
  hPhoPt->Scale(1./scale);
  hMllg->Scale(1./scale);

  cout << "number of pho: " << ipho << endl;

  file.Close();

  fout->cd();
  hPhoPt->Write();
  hMllg->Write();
  htotwei->Write();

  outtree->Write();

  cout << "write to file" << endl;
  fout->Write();
  fout->Close();

  cout << "done !!!!!!!!!" << endl;
  
}
