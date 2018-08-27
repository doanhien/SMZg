#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"


void correlationPlot() {

  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.3, "Y");

  TFile *fin = new TFile("../ana/minitrees_Zee/DYJets_amcatnlo_TnP_Zee_chisocor_NoGenMatch_PhoPresel_Corr_BDTUpto6000.root", "read");

  TTree *t = (TTree*) fin->Get("passingIdTree");

  TProfile *hsieie_sieip = new TProfile("hsieie_sieip","sigma ietaieta vs ietaiphi", 50, 0., 0.05, -0.001, 0.001);
  TProfile *hsieie_r9 = new TProfile("hsieie_r9","sigma ietaieta vs R9", 50, 0., 0.05, 0, 1.);
  TProfile *hsieie_s4 = new TProfile("hsieie_s4","sigma ietaieta vs S4", 50, 0., 0.05, 0, 1.);
  TProfile *hsieie_etaw = new TProfile("hsieie_etaw","sigma ietaieta vs SCEta width", 50, 0., 0.05, 0, 1.);
  TProfile *hsieie_phiw = new TProfile("hsieie_phiw","sigma ietaieta vs SCEta width", 50, 0., 0.05, 0, 1.);
  TProfile *hsieie_rho = new TProfile("hsieie_rho","sigma ietaieta vs SCEta width", 50, 0., 0.05, 0, 50);

  TProfile *hsieip_r9 = new TProfile("hsieip_r9","sigma ietaieta vs R9", 50, -1e-3, 1e-3, 0, 1.);
  TProfile *hsieip_s4 = new TProfile("hsieip_s4","sigma ietaieta vs S4", 50, -1e-3, 1e-3, 0, 1.);
  TProfile *hsieip_etaw = new TProfile("hsieip_etaw","sigma ietaieta vs SCEta width", 50, -1e-3, 1e-3, 0, 1.);
  TProfile *hsieip_phiw = new TProfile("hsieip_phiw","sigma ietaieta vs SCEta width", 50, -1e-3, 1e-3, 0, 1.);
  TProfile *hsieip_rho = new TProfile("hsieip_rho","sigma ietaieta vs SCEta width", 50, -1e-3, 1e-3, 0, 50);

  TProfile *hr9_s4 = new TProfile("hr9_s4","sigma ietaieta vs S4", 50, 0, 1., 0, 1.);
  TProfile *hr9_etaw = new TProfile("hr9_etaw","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 1.);
  TProfile *hr9_phiw = new TProfile("hr9_phiw","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 1.);
  TProfile *hr9_rho = new TProfile("hr9_rho","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 50);

  TProfile *hs4_etaw = new TProfile("hs4_etaw","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 1.);
  TProfile *hs4_phiw = new TProfile("hs4_phiw","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 1.);
  TProfile *hs4_rho = new TProfile("hs4_rho","sigma ietaieta vs SCEta width", 50, 0, 1., 0, 50);

  TProfile *hetaw_phiw = new TProfile("hetaw_phiw","sigma ietaieta vs SCEta width", 50, 0, 0.2, 0, 1.);
  TProfile *hetaw_rho = new TProfile("hetaw_rho","sigma ietaieta vs SCEta width", 50, 0, 0.2, 0, 50);


  TCut cut = "Zm>70 && Zm<110 && Probe_Pt>15 && passPhoID_Zg==1";

  t->Draw("Probe_sieip:Probe_sieie >> hsieie_sieip", cut);
  t->Draw("Probe_R9:Probe_sieie >> hsieie_r9", cut);
  t->Draw("Probe_s4Full5x5:Probe_sieie >> hsieie_s4", cut);
  t->Draw("Probe_SCEtaWidth_rw:Probe_sieie_rw >> hsieie_etaw", cut);
  t->Draw("Probe_scphiwidth:Probe_sieie >> hsieie_phiw", cut);
  t->Draw("rho:Probe_sieie >> hsieie_rho", cut);

  t->Draw("Probe_R9:Probe_sieip >> hsieip_r9", cut);
  t->Draw("Probe_s4Full5x5:Probe_sieip >> hsieip_s4", cut);
  t->Draw("Probe_SCEtaWidth:Probe_sieip >> hsieip_etaw", cut);
  t->Draw("Probe_scphiwidth:Probe_sieip >> hsieip_phiw", cut);
  t->Draw("rho:Probe_sieip >> hsieip_rho", cut);

  t->Draw("Probe_s4Full5x5:Probe_R9 >> hr9_s4", cut);
  t->Draw("Probe_SCEtaWidth:Probe_R9 >> hr9_etaw", cut);
  t->Draw("Probe_scphiwidth:Probe_R9 >> hr9_phiw", cut);
  t->Draw("rho:Probe_R9 >> hr9_rho", cut);

  t->Draw("Probe_SCEtaWidth:Probe_s4Full5x5 >> hs4_etaw", cut);
  t->Draw("Probe_scphiwidth:Probe_s4Full5x5 >> hs4_phiw", cut);
  t->Draw("rho:Probe_s4Full5x5 >> hs4_rho", cut);

  t->Draw("Probe_scphiwidth:Probe_SCEtaWidth >> hetaw_phiw", cut);
  t->Draw("rho:Probe_SCEtaWidth >> hetaw_rho", cut);


  TCanvas *c1 = new TCanvas("c1", "c1", 1050, 800);
  c1->Divide(2,2);
  c1->cd(1);
  hsieie_sieip->GetYaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieie_sieip->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_sieip->Draw();

  c1->cd(2);
  hsieie_r9->GetYaxis()->SetTitle("R9");
  hsieie_r9->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_r9->Draw();

  c1->cd(3);
  hsieie_s4->GetYaxis()->SetTitle("E2x2/E5x5");
  hsieie_s4->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_s4->Draw();

  c1->cd(4);
  hsieie_etaw->GetYaxis()->SetTitle("SCEta width");
  hsieie_etaw->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_etaw->Draw();

  TCanvas *c2 = new TCanvas("c2", "c2", 1050, 400);
  c2->Divide(2,1);
  c2->cd(1);
  hsieie_phiw->GetYaxis()->SetTitle("SCPhi width");
  hsieie_phiw->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_phiw->Draw();

  c2->cd(2);
  hsieie_rho->GetYaxis()->SetTitle("#rho");
  hsieie_rho->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_rho->Draw();


  //for sigma_ieta_ipho
  TCanvas *c3 = new TCanvas("c3", "c3", 1050, 800);
  c3->Divide(2,2);
  c3->cd(1);
  hsieip_r9->GetXaxis()->SetNdivisions(505);
  hsieip_r9->GetYaxis()->SetTitle("R9");
  hsieip_r9->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_r9->Draw();

  c3->cd(2);
  hsieip_s4->GetXaxis()->SetNdivisions(505);
  hsieip_s4->GetYaxis()->SetTitle("E2x2/E5x5");
  hsieip_s4->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_s4->Draw();

  c3->cd(3);
  hsieip_etaw->GetXaxis()->SetNdivisions(505);
  hsieip_etaw->GetYaxis()->SetTitle("SCEta width");
  hsieip_etaw->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_etaw->Draw();

  c3->cd(4);
  hsieip_phiw->GetXaxis()->SetNdivisions(505);
  hsieip_phiw->GetYaxis()->SetTitle("SCPhi width");
  hsieip_phiw->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_phiw->Draw();


  TCanvas *c4 = new TCanvas("c4", "c4", 550, 400);
  c4->cd();
  hsieip_rho->GetXaxis()->SetNdivisions(505);
  hsieip_rho->GetYaxis()->SetTitle("#rho");
  hsieip_rho->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_rho->Draw();

  //for r9
  TCanvas *cr9 = new TCanvas("cr9", "cr9", 1050, 800);
  cr9->Divide(2,2);
  cr9->cd(1);
  hr9_s4->GetYaxis()->SetTitle("E2x2/E5x5");
  hr9_s4->GetXaxis()->SetTitle("R9");
  hr9_s4->Draw();

  cr9->cd(2);
  hr9_etaw->GetYaxis()->SetTitle("SCEta width");
  hr9_etaw->GetXaxis()->SetTitle("R9");
  hr9_etaw->Draw();

  cr9->cd(3);
  hr9_phiw->GetYaxis()->SetTitle("SCPhi width");
  hr9_phiw->GetXaxis()->SetTitle("R9");
  hr9_phiw->Draw();

  cr9->cd(4);
  hr9_rho->GetYaxis()->SetTitle("rho");
  hr9_rho->GetXaxis()->SetTitle("R9");
  hr9_rho->Draw();

  //for s4
  TCanvas *cs4 = new TCanvas("cs4", "cs4", 1050, 800);
  cs4->Divide(2,2);
  cs4->cd(1);
  hs4_etaw->GetYaxis()->SetTitle("SCEta width");
  hs4_etaw->GetXaxis()->SetTitle("E2x2/E5x5");
  hs4_etaw->Draw();

  cs4->cd(2);
  hs4_phiw->GetYaxis()->SetTitle("SCPhi width");
  hs4_phiw->GetXaxis()->SetTitle("E2x2/E5x5");
  hs4_phiw->Draw();

  cs4->cd(3);
  hs4_rho->GetYaxis()->SetTitle("rho");
  hs4_rho->GetXaxis()->SetTitle("E2x2/E5x5");
  hs4_rho->Draw();


  //for etaw
  TCanvas *cetaw = new TCanvas("cetaw", "cetaw", 1050, 400);
  cetaw->Divide(2,1);
  cetaw->cd(1);
  hetaw_phiw->GetYaxis()->SetTitle("SCPhi width");
  hetaw_phiw->GetXaxis()->SetTitle("SCEta width");
  hetaw_phiw->Draw();

  cetaw->cd(2);
  hetaw_rho->GetYaxis()->SetTitle("#rho");
  hetaw_rho->GetXaxis()->SetTitle("SCEta width");
  hetaw_rho->Draw();

  /*
  c1->SaveAs("plots/correlation_InputVarTrain/sieie_1.pdf");
  c2->SaveAs("plots/correlation_InputVarTrain/sieie_2.pdf");
  c3->SaveAs("plots/correlation_InputVarTrain/sieip_1.pdf");
  c4->SaveAs("plots/correlation_InputVarTrain/sieip_2.pdf");
  cr9->SaveAs("plots/correlation_InputVarTrain/r9.pdf");
  cs4->SaveAs("plots/correlation_InputVarTrain/s4.pdf");
  cetaw->SaveAs("plots/correlation_InputVarTrain/etawidth.pdf");
  */

}
