#include <TFile.h>
#include <TH1.h>
#include <TStyle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>
#include <TCut.h>

#include <iostream>

using namespace std;

void histStyle(TH1F *h1) {
  h1->SetTitleSize(0.14, "XYZ");
  h1->SetLabelSize(0.14, "XYZ");
  h1->SetTitleFont(42, "XYZ");
  h1->GetYaxis()->SetTitle("data/ MC");
  //h1->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");
  h1->GetYaxis()->SetTitleOffset(0.55);
  h1->GetXaxis()->SetTitleOffset(1.1);
  h1->GetYaxis()->SetRangeUser(0., 2);
  h1->GetYaxis()->SetNdivisions(5);
  //h1->Draw("ep");
}

void histStyle1 (TH1F *h) {
  h->SetFillColor(kCyan-7);
  h->SetLineColor(kCyan-7);
}

void histStyle2 (TH1F *h) {
  h->SetFillColor(kBlue);
  h->SetLineColor(kBlue);
}

void histStyle3 (TH1F *h) {
  h->SetFillColor(kBlue-9);
  h->SetLineColor(kBlue-9);
}

void histStyle4 (TH1F *h) {
  h->SetFillColor(kViolet+1);
  h->SetLineColor(kViolet+1);

}


void IdPlot(bool ele = true, TString hname = "hlepd0_lead") {

  gStyle->SetOptStat(0);

  cout << "get started" << endl;
  TString infname;
  if (ele) infname  = "histo/Id/Ele_EB_";
  else infname  = "histo/Id/Mu_EB_";

  const int nFile = 11;
  TFile *fin[nFile];
  TH1F *h[nFile];

  TH1F  *hMC = new TH1F();

  if (ele)
    fin[0] = new TFile(infname + "DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_NoPhoSel.root", "read");
  else fin[0] = new TFile(infname + "DoubleMu_Run2016_FebReminiAOD_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[1] = new TFile(infname + "Zg_aMCatNLO_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[2] = new TFile(infname + "ZJets_aMCatNLO_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[3] = new TFile(infname + "TT_Powheg_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[4] = new TFile(infname + "WWTo2L2Nu_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[5] = new TFile(infname + "WWToLNuQQ_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[6] = new TFile(infname + "WZTo3LNu_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[7] = new TFile(infname + "WZTo2L2Q_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[8] = new TFile(infname + "ZZTo2L2Nu_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[9] = new TFile(infname + "ZZTo2L2Q_Summer16_TMVA420_NoPhoSel.root", "read");
  fin[10] = new TFile(infname + "ZZTo4L_Summer16_TMVA420_NoPhoSel.root", "read");


  cout << "done opening files" << endl;
  for (int i = 0; i < nFile; i++) {
    h[i] = (TH1F*) fin[i]->Get(hname.Data());
  }

  for (int i = 1; i < nFile; i++) {
    if (i == 1) hMC = (TH1F*) h[i]->Clone();
    else hMC->Add(h[i]);
  }

  cout << "total entries from mc: " << hMC->Integral() << endl;
  cout << "total entries from data: " << h[0]->Integral() << endl;
  cout << "total entries from Zg: " << h[1]->Integral() << endl;
  cout << "total entries from Zj: " << h[2]->Integral() << endl;
  cout << "total entries from tt: " << h[3]->Integral() << endl;

  //htest = (TH1F*) fin[3]->Get(hname.Data());
  //cout << "testing event: " << htest->Integral() << endl;
  h[4]->Add(h[5]);
  h[4]->Add(h[6]);
  h[4]->Add(h[7]);
  h[4]->Add(h[8]);
  h[4]->Add(h[9]);
  h[4]->Add(h[10]);

  cout << "total entries from VV: " << h[4]->Integral() << endl;

  h[0]->Sumw2();
  TH1F *hRatio = (TH1F*) h[0]->Clone();
  hRatio->Divide(hMC);
  histStyle(hRatio);

  histStyle1(h[1]);
  histStyle2(h[2]);
  histStyle3(h[3]);
  histStyle4(h[4]);

  THStack *hstack = new THStack("hstack", "");
  hstack->Add(h[4]);
  hstack->Add(h[3]);
  hstack->Add(h[2]);
  hstack->Add(h[1]);

  TString xname, yname, outname; 

  if (ele) {
    if (hname == "hlepd0_lead") {
      xname = "D0";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_d0";}
    if (hname == "hlepdz_lead") {
      xname = "DZ";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_dz";}
    if (hname == "hSIP_lead") {
      xname = "SIP";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_SIP"; }
    if (hname == "hchiso_lead") {
      xname = "charged iso";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_chargedIso";  }
    if (hname == "hneuiso_lead") {
      xname = "neutral iso";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_neuIso"; }
    if (hname == "hphoiso_lead") {
      xname = "photon iso";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_phoIso"; }
    if (hname == "hmva_lead") {
      xname = "MVA";   yname = "Events"; outname = "plots/Id/data_MC_leadEle_mva"; }

    if (hname == "hlepd0_trail") {
      xname = "D0";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_d0"; }
    if (hname == "hlepdz_trail") {
      xname = "DZ";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_dz"; }
    if (hname == "hSIP_trail") {
      xname = "SIP";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_SIP"; }
    if (hname == "hchiso_trail") {
      xname = "charged iso";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_chargedIso"; }
    if (hname == "hneuiso_trail") {
      xname = "neutral iso";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_neuIso"; }
    if (hname == "hphoiso_trail") {
      xname = "photon iso";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_phoIso"; }
    if (hname == "hmva_trail") {
      xname = "MVA";   yname = "Events"; outname = "plots/Id/data_MC_trailEle_mva"; }
    if (hname == "hphoHoE") {
      xname = "photon H/E";   yname = "Events"; outname = "plots/Id/data_MC_Ele_phoHoE"; }
    if (hname == "hphosieie") {
      xname = "photon #sigma_{i#etai#eta}";   yname = "Events"; outname = "plots/Id/data_MC_Ele_phoSieie"; }
    if (hname == "hphochiso") {
      xname = "charged iso^{#gamma}";   yname = "Events"; outname = "plots/Id/data_MC_Ele_phoChargedIso"; }
    else if (hname == "hphowchiso") {
      xname = "worst charged iso^{#gamma}";   yname = "Events"; outname = "plots/Id/data_MC_Ele_phoWorstChargedIso"; }
    if (hname == "hsieip") {
      xname = "photon #sigma_{i#etai#phi}"; yname = "Events"; outname = "plots/Id/data_MC_Ele_phoSieip";}
    if (hname == "hetawidth") {
      xname = "#eta_{SC}^{width}"; yname = "Events"; outname ="plots/Id/data_MC_Ele_SCEtaWidth";}
    if (hname == "hphiwidth") {
      xname = "#phi_{SC}^{width}"; yname = "Events"; outname ="plots/Id/data_MC_Ele_SCPhiWidth";}
    if (hname == "hscRawE") {
      xname = "E^{raw}_{SC} [GeV]"; yname = "Events"; outname ="plots/Id/data_MC_Ele_SCRawE";}
    if (hname == "hs4") {
      xname = "E_{2x2}/E_{5x5}"; yname = "Events"; outname ="plots/Id/data_MC_Ele_S4";}
    if (hname == "hr9") {
      xname = "photon R9"; yname = "Events"; outname ="plots/Id/data_MC_Ele_phoR9";}
    if (hname == "hphophi") {
      xname = "#phi^{#gamma}"; yname = "Events"; outname ="plots/Id/data_MC_Ele_PhoPhi";}
    if (hname == "hrho") {
      xname = "#rho"; yname = "Events"; outname ="plots/Id/data_MC_Ele_rho";}
    if (hname == "hphossmva") {
      xname = "showershape MVA";   yname = "Events/0.1"; outname += "plots/Id/Ele_ssmva"; }


  }
  else {
    if (hname == "hlepd0_lead") {
      xname = "D0";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_d0"; }
    if (hname == "hlepdz_lead") {
      xname = "DZ";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_dz"; }
    if (hname == "hSIP_lead") {
      xname = "SIP";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_SIP"; }
    if (hname == "hchiso_lead") {
      xname = "charged iso";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_chargedIso"; }
    if (hname == "hneuiso_lead") {
      xname = "neutral iso";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_neuIso"; }
    if (hname == "hphoiso_lead") {
      xname = "photon iso";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_phoIso"; }
    if (hname == "hsigEOverE_lead") {
      xname = "#sigma_{pT}^{trk}/p_{T}^{trk}";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_sigPtOverPt"; }
    if (hname == "hmuSta") {
      xname = "muon stations";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_station"; }
    if (hname == "hmuPixhit") {
      xname = "muon pixel hits";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_PixHits"; }
    if (hname == "hmutrkLayer") {
      xname = "muon tracker layers";   yname = "Events"; outname = "plots/Id/data_MC_leadMu_trkLayer"; }

    if (hname == "hlepd0_trail") {
      xname = "D0";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_d0"; }
    if (hname == "hlepdz_trail") {
      xname = "DZ";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_dz"; }
    if (hname == "hSIP_trail") {
      xname = "SIP";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_SIP"; }
    if (hname == "hchiso_trail") {
      xname = "charged iso";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_chargedIso"; }
    if (hname == "hneuiso_trail") {
      xname = "neutral iso";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_neuIso"; }
    if (hname == "hphoiso_trail") {
      xname = "photon iso";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_phoIso"; }
    if (hname == "hsigEOverE_trail")  {
      xname = "#sigma_{pT}^{trk}/p_{T}^{trk}";   yname = "Events"; outname = "plots/Id/data_MC_trailMu_sigPtOverPt"; }

    if (hname == "hphoHoE") {
      xname = "photon H/E";   yname = "Events"; outname = "plots/Id/data_MC_Mu_phoHoE"; }
    if (hname == "hphosieie") {
      xname = "photon #sigma_{i#etai#eta}";   yname = "Events"; outname = "plots/Id/data_MC_Mu_phoSieie"; }
    if (hname == "hphochiso") {
      xname = "charged iso^{#gamma}";   yname = "Events"; outname = "plots/Id/data_MC_Mu_phoChargedIso"; }
    if (hname == "hphowchiso") {
      xname = "worst charged iso^{#gamma}";   yname = "Events"; outname = "plots/Id/data_MC_Mu_phoWorstChargedIso"; }
    if (hname == "hsieip") {
      xname = "photon #sigma_{i#etai#phi}"; yname = "Events"; outname = "plots/Id/data_MC_Mu_phoSieip";}
    if (hname == "hetawidth") {
      xname = "#eta_{SC}^{width}"; yname = "Events"; outname ="plots/Id/data_MC_Mu_SCEtaWidth";}
    if (hname == "hphiwidth") {
      xname = "#phi_{SC}^{width}"; yname = "Events"; outname ="plots/Id/data_MC_Mu_SCPhiWidth";}
    if (hname == "hscRawE") {
      xname = "E^{raw}_{SC} [GeV]"; yname = "Events"; outname ="plots/Id/data_MC_Mu_SCRawE";}
    if (hname == "hs4") {
      xname = "E_{2x2}/E_{5x5}"; yname = "Events"; outname ="plots/Id/data_MC_Mu_S4";}
    if (hname == "hr9") {
      xname = "photon R9"; yname = "Events"; outname ="plots/Id/data_MC_Mu_phoR9";}
    if (hname == "hphophi") {
      xname = "#phi^{#gamma}"; yname = "Events"; outname ="plots/Id/data_MC_Mu_PhoPhi";}
    if (hname == "hrho") {
      xname = "#rho"; yname = "Events"; outname ="plots/Id/data_MC_Mu_rho";}
    if (hname == "hphossmva") {
      xname = "showershape MVA";   yname = "Events/0.1"; outname += "plots/Id/Ele_ssmva"; }

  }

  //cout << "outname for saving: " << outname.Data() << endl;
  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.05);
  tx.SetTextFont(52);

  //plotting for canvas
  TCanvas* c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  c1->SetBottomMargin(0.15);
  c1->SetLogy();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  if (hname == "hphoHoE") pad1->SetLogy();
  h[0]->GetYaxis()->SetRangeUser(0.9, h[0]->GetMaximum()*2);
  h[0]->GetXaxis()->SetLabelSize(0);
  h[0]->GetYaxis()->SetTitle(yname);
  h[0]->Draw();
  hstack->Draw("histsame");
  h[0]->Draw("same");
  pad1->RedrawAxis();
  TLegend *lg = new TLegend(0.68, 0.68, 0.85, 0.88);
  lg->AddEntry(h[0], "data", "pe");
  lg->AddEntry(h[1], "Z#gamma", "f");
  lg->AddEntry(h[2], "DYJets", "f");
  lg->AddEntry(h[3], "T#barT", "f");
  lg->AddEntry(h[4], "VV", "f");
  lg->SetLineColor(0);
  lg->SetFillColor(0);
  lg->SetShadowColor(0);
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->Draw();
  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  hRatio->GetXaxis()->SetTickLength(0.05);
  hRatio->GetXaxis()->SetTitle(xname);
  hRatio->Draw("ep");

  c1->SaveAs(outname + "_EB_NoPhoSel.pdf");

}
