#include "TFile.h"
#include "TTree.h"
#include "TH2.h"

void draw(bool barrel = true) {

  gStyle->SetOptStat(0);

  /*

// ---- plot mva vs charged iso ----------//
  TFile *fsig = new TFile("../ana/minitrees/Zg_aMCatNLO_Summer16_TMVA420_LooseSieie_chisocor.root", "read");
  TFile *fbkg = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_LooseSieie_chisocor_genmatch.root", "read");


  TTree *tsig = (TTree*) fsig->Get("outtree");
  TTree *tbkg = (TTree*) fbkg->Get("outtree");

  TH2F *hmva_iso_sig = new TH2F("hmva_iso_sig", "mva vs iso", 20, -1, 1, 15, 0, 15);
  TH2F *hmva_iso_bkg = new TH2F("hmva_iso_bkg", "mva vs iso", 20, -1, 1, 15, 0, 15);

  TCut cut = "z_charge==0 && trig_Ele23_Ele12==1 && leptType==11 && lept1_pt>20 && fabs(lept0_eta)<2.5 && fabs(lept1_eta)<2.5";
  //TCut cut = "z_charge==0 && trig_Mu17_Mu8==1 && leptType==13 && lept1_pt>20";

  tsig->Draw("gamma_ChIso:gamma_ssmva >> hmva_iso_sig", cut);
  tbkg->Draw("gamma_ChIso:gamma_ssmva >> hmva_iso_bkg", cut);

  cs = new TCanvas("cs", "cs", 680, 600);
  cs->cd();
  cs->SetRightMargin(0.17);
  cs->SetLeftMargin(0.14);
  hmva_iso_sig->GetYaxis()->SetTitle("charged isolation [GeV]");
  hmva_iso_sig->GetXaxis()->SetTitle("photon MVA");
  hmva_iso_sig->Draw("colz");

  cb = new TCanvas("cb", "cb", 680, 600);
  cb->cd();
  cb->SetRightMargin(0.17);
  cb->SetLeftMargin(0.14);
  hmva_iso_bkg->GetYaxis()->SetTitle("charged isolation [GeV]");
  hmva_iso_bkg->GetXaxis()->SetTitle("photon MVA");
  hmva_iso_bkg->Draw("colz");

  cs->SaveAs("plots/Id/mva_vs_chargedIso_signalZg_Ele.pdf");
  cb->SaveAs("plots/Id/mva_vs_chargedIso_bkgZJet_Ele.pdf");
  */


  // ----- mva between MC and data for Zee---------//

  TFile *fda = new TFile("../ana/minitrees_Zee/SingleEle_Run2016_TnP_Zee_chisocor_NewSSCor.root", "read");
  TFile *fDY = new  TFile("../ana/minitrees_Zee/DYJets_amcatnlo_TnP_Zee_chisocor_NewSSCor.root", "read");
  TFile *fttb = new TFile("../ana/minitrees_Zee/TTbar_TnP_Zee_chisocor_NewSSCor.root", "read");
  TFile *fwj = new TFile("../ana/minitrees_Zee/WJetsToLNu_TnP_Zee_chisocor_NewSSCor.root", "read");

  TTree *tda = (TTree*) fda->Get("passingIdTree");
  TTree *tDY = (TTree*) fDY->Get("passingIdTree");
  TTree *tttb = (TTree*) fttb->Get("passingIdTree");
  TTree *twj = (TTree*) fwj->Get("passingIdTree");

  TH1F *hmva_da = new TH1F("hmva_da", "", 20, -1, 1);
  TH1F *hmva_DY =new TH1F("hmva_DY", "", 20, -1, 1);
  TH1F *hmva_ttb =new TH1F("hmva_ttb", "", 20, -1, 1);
  TH1F *hmva_wj =new TH1F("hmva_wj", "", 20, -1, 1);

  TCut cut = "Zm>70 && Zm<110 && Probe_Pt>15 && passPhoID_Zg==1";

  //if (barrel) cut += "fabs(Probe_SCEta)<1.4442 && Probe_PhoChIso<2.";
  //else cut += "fabs(Probe_SCEta)>1.566 && Probe_PhoChIso<1.5";
  if (barrel) cut += "fabs(Probe_SCEta)<1.5";
  else cut += "fabs(Probe_SCEta)>1.5";

  TCut mcwei = "puweigj_65nb*genWeight*Tag_RecoSF*Tag_SelSF*Probe_phoSF";
  //TCut mcwei = "puweigj_65nb*genWeight*Tag_RecoSF*Tag_SelSF";


  tda->Draw("Probe_ssmva >> hmva_da", cut);
  tDY->Draw("Probe_ssmva >> hmva_DY", cut*mcwei);
  tttb->Draw("Probe_ssmva >> hmva_ttb", cut*mcwei);
  twj->Draw("Probe_ssmva >> hmva_wj", cut*mcwei);

  float w_dy = 2.60897;
  float w_tt = 0.0399477;
  float w_wj = 35900*6.072e+04/160712265; //effective event = 160712265

  hmva_DY->Scale(w_dy);
  hmva_ttb->Scale(w_tt);
  hmva_wj->Scale(w_wj);

  cout << "Data: "  << hmva_da->Integral() << endl;
  cout << "DY: "    << hmva_DY->Integral() << endl;
  cout << "TTbar: " << hmva_ttb->Integral()<< endl;
  cout << "WJet: "  << hmva_wj->Integral() << endl;

  hmva_da->Sumw2();
  hmva_DY->SetLineColor(kOrange-2);
  hmva_DY->SetFillColor(kOrange-2);
  hmva_ttb->SetLineColor(kRed-6);
  hmva_ttb->SetFillColor(kRed-6);
  hmva_wj->SetLineColor(kCyan-3);
  hmva_wj->SetFillColor(kCyan-3);

  THStack *hmc = new THStack("hmc", "");
  hmc->Add(hmva_wj);
  hmc->Add(hmva_ttb);
  hmc->Add(hmva_DY);

  //hmva_da->Scale(1./hmva_da->Integral());
  //hmva_DY->Scale(1./hmva_DY->Integral());

  c1 = new TCanvas("c1", "plot of mva", 650, 650);
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  //pad1->SetLogy();
  pad1->SetFillColor(0);

  float ymax = hmc->GetMaximum();
  if(hmva_da->GetMaximum() > ymax) ymax = hmva_da->GetMaximum();
  hmva_da->SetMaximum(ymax*1.2);

  hmva_da->GetYaxis()->SetTitleOffset(1.6);
  hmva_da->SetTitleSize(0.05, "XYZ");
  hmva_da->SetLabelSize(0.05, "XYZ");
  hmva_da->SetTitleFont(42, "XYZ");
  hmva_da->GetYaxis()->SetTitle("Entries/0.1");

  hmva_da->SetLabelSize(0, "X");
  hmva_da->SetTitleOffset(1.5, "Y");
  hmva_da->SetTitleOffset(1.3, "X");
  //hmva_da->GetYaxis()->SetRangeUser(0.1, hmva_da->GetMaximum()*10);
  hmva_da->Draw();
  //hmva_DY->Draw("hist same");
  hmc->Draw("histsame");
  hmva_da->Draw("same");
  TLegend *lg = new TLegend(0.28, 0.5, 0.45, 0.78);
  lg->AddEntry(hmva_da, "data", "pe");
  lg->AddEntry(hmva_DY, "DY","f");
  lg->AddEntry(hmva_ttb, "TT#bar","f");
  lg->AddEntry(hmva_wj, "WJets","f");
  lg->Draw();
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  pad1->RedrawAxis();


  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.04);
  tx.SetTextFont(42);
  tx.DrawLatex(0.42, 0.84, "Z#gamma#rightarrow ee#gamma");
  if (barrel) tx.DrawLatex(0.4, 0.78, "0 < |#eta^{#gamma}| < 1.5");
  else tx.DrawLatex(0.4,0.78, "1.5 < |#eta^{#gamma}| < 2.5");


  c1->cd();
  TH1F *hratio = (TH1F*) hmva_da->Clone();
  hmva_DY->Add(hmva_ttb);
  hmva_DY->Add(hmva_wj);
  hratio->Divide(hmva_DY);
  hratio->SetLineColor(1);

  cout << "DY: "    << hmva_DY->Integral() << endl;

  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();

  hratio->GetYaxis()->SetTitleOffset(1.6);
  hratio->SetMinimum(-6);
  hratio->SetTitleSize(0.14, "XYZ");
  hratio->SetLabelSize(0.14, "XYZ");
  hratio->SetTitleFont(42, "XYZ");
  hratio->GetYaxis()->SetTitle("data/ MC");
  hratio->GetXaxis()->SetTitle("photon MVA");
  hratio->GetYaxis()->SetTitleOffset(0.55);
  hratio->GetXaxis()->SetTitleOffset(1.1);
  hratio->GetYaxis()->SetRangeUser(0., 2.);
  hratio->GetYaxis()->SetNdivisions(8);
  hratio->Draw();

  TLine *l = new TLine(-1,1,1,1);
  l->SetLineWidth(2);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");

  //if (barrel)
  //c1->SaveAs("plots/Id/Zee/photonMVA_Zee_EB_PhoPreSel.pdf");
  //else c1->SaveAs("plots/Id/Zee/photonMVA_Zee_EE_PhoPreSel.pdf");



}
