#include "TH1.h"

void compare_private_combine() {

  //private code
  float yield_eb_ele_pri[] = {9291.2, 5066.86, 2671.98, 1484.16, 1535.51};
  float yield_ee_ele_pri[] = {3460.69, 1620.48, 832.388, 525.776, 572.369};
  float yield_err_eb_ele_pri[] = {131.014, 89.675, 63.4491, 46.7399, 47.6311};
  float yield_err_ee_ele_pri[] = {81.4716, 52.9539, 37.192, 27.6333, 28.2478};

  float yield_eb_mu_pri[] = {20926.3, 11648.1, 6207.89, 3003.91, 3086.19};
  float yield_ee_mu_pri[] = {8081.41, 3947.44, 2070.19, 1033.03, 1149.43};
  float yield_err_eb_mu_pri[] = {195.153, 134.485, 95.4255, 67.4658, 67.6339};
  float yield_err_ee_mu_pri[] = {124.46, 81.0706, 56.5616, 39.34, 40.199};

  //combine tool
  float yield_eb_ele_comb[] = {9291.43, 5066.87, 2672.03, 1484.22, 1535.54};
  float yield_ee_ele_comb[] = {3460.67, 1620.48, 832.387, 525.787, 572.367};
  float yield_err_eb_ele_comb[] = {136.245, 93.7407, 66.4127, 48.9583, 49.8353};
  float yield_err_ee_ele_comb[] = {85.1231, 55.4341, 38.9796, 29.1184, 29.7928};

  float yield_eb_mu_comb[] = {20926.4, 11648.1, 6207.72, 3003.9, 3086.18};
  float yield_ee_mu_comb[] = {8081.41, 3947.5, 2070.27, 1033.04, 1149.45};
  float yield_err_eb_mu_comb[] = {203.057, 140.588, 99.8901, 70.5602, 70.7941};
  float yield_err_ee_mu_comb[] = {129.768, 84.8848, 59.5033, 41.3974, 42.4039};

  int nbin = 5;
  float ptbin[] = {15, 20, 25, 30, 35, 45};
  TH1F *hyield_eb_ele_pri = new TH1F("hyield_eb_ele_pri", "yield of EB ele private", 5, ptbin);
  TH1F *hyield_ee_ele_pri = new TH1F("hyield_ee_ele_pri", "yield of EB ele private", 5, ptbin);
  TH1F *hyield_eb_mu_pri  = new TH1F("hyield_eb_mu_pri", "yield of EB ele private", 5, ptbin);
  TH1F *hyield_ee_mu_pri  = new TH1F("hyield_ee_mu_pri", "yield of EB ele private", 5, ptbin);

  TH1F *hyield_eb_ele_comb = new TH1F("hyield_eb_ele_comb", "yield of EB ele combine", 5, ptbin);
  TH1F *hyield_ee_ele_comb = new TH1F("hyield_ee_ele_comb", "yield of EB ele combine", 5, ptbin);
  TH1F *hyield_eb_mu_comb  = new TH1F("hyield_eb_mu_comb", "yield of EB ele combine", 5, ptbin);
  TH1F *hyield_ee_mu_comb  = new TH1F("hyield_ee_mu_comb", "yield of EB ele combine", 5, ptbin);

  for (int i = 1; i <= nbin; i++) {
    
    //cout << "processing " << i << " bin" << endl;
    hyield_eb_ele_pri->SetBinContent(i, yield_eb_ele_pri[i-1]);
    hyield_ee_ele_pri->SetBinContent(i, yield_ee_ele_pri[i-1]);
    hyield_eb_mu_pri->SetBinContent(i, yield_eb_mu_pri[i-1]);
    hyield_ee_mu_pri->SetBinContent(i, yield_ee_mu_pri[i-1]);

    hyield_eb_ele_pri->SetBinError(i, yield_err_eb_ele_pri[i-1]);
    hyield_ee_ele_pri->SetBinError(i, yield_err_ee_ele_pri[i-1]);
    hyield_eb_mu_pri->SetBinError(i, yield_err_eb_mu_pri[i-1]);
    hyield_ee_mu_pri->SetBinError(i, yield_err_ee_mu_pri[i-1]);


    hyield_eb_ele_comb->SetBinContent(i, yield_eb_ele_comb[i-1]);
    hyield_ee_ele_comb->SetBinContent(i, yield_ee_ele_comb[i-1]);
    hyield_eb_mu_comb->SetBinContent(i, yield_eb_mu_comb[i-1]);
    hyield_ee_mu_comb->SetBinContent(i, yield_ee_mu_comb[i-1]);

    hyield_eb_ele_comb->SetBinError(i, yield_err_eb_ele_comb[i-1]);
    hyield_ee_ele_comb->SetBinError(i, yield_err_ee_ele_comb[i-1]);
    hyield_eb_mu_comb->SetBinError(i, yield_err_eb_mu_comb[i-1]);
    hyield_ee_mu_comb->SetBinError(i, yield_err_ee_mu_comb[i-1]);
  }


  hyield_eb_ele_pri->SetMarkerStyle(20);
  hyield_eb_ele_pri->SetMarkerColor(2);
  hyield_eb_ele_pri->SetLineColor(2);

  hyield_ee_ele_pri->SetMarkerStyle(20);
  hyield_ee_ele_pri->SetMarkerColor(2);
  hyield_ee_ele_pri->SetLineColor(2);

  hyield_eb_mu_pri->SetMarkerStyle(20);
  hyield_eb_mu_pri->SetMarkerColor(2);
  hyield_eb_mu_pri->SetLineColor(2);

  hyield_ee_mu_pri->SetMarkerStyle(20);
  hyield_ee_mu_pri->SetMarkerColor(2);
  hyield_ee_mu_pri->SetLineColor(2);

  hyield_eb_ele_comb->SetMarkerStyle(24);
  hyield_eb_ele_comb->SetMarkerColor(4);
  hyield_eb_ele_comb->SetLineColor(4);

  hyield_ee_ele_comb->SetMarkerStyle(24);
  hyield_ee_ele_comb->SetMarkerColor(4);
  hyield_ee_ele_comb->SetLineColor(4);

  hyield_eb_mu_comb->SetMarkerStyle(24);
  hyield_eb_mu_comb->SetMarkerColor(4);
  hyield_eb_mu_comb->SetLineColor(4);

  hyield_ee_mu_comb->SetMarkerStyle(24);
  hyield_ee_mu_comb->SetMarkerColor(4);
  hyield_ee_mu_comb->SetLineColor(4);

  //hyield_eb_ele_comb->Print("all");

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
  c1->cd();
  TPad *pad11 = new TPad("pad11", "pad11", 0, 0.25, 1, 1);
  pad11->SetBottomMargin(0.02);
  pad11->Draw();
  pad11->cd();
  //pad11->SetLogy();
  pad11->SetFillColor(0);

  hyield_eb_mu_comb->GetYaxis()->SetTitleOffset(1.6);
  hyield_eb_mu_comb->SetTitleSize(0.05, "XYZ");
  hyield_eb_mu_comb->SetLabelSize(0.05, "XYZ");
  hyield_eb_mu_comb->SetTitleFont(42, "XYZ");
  hyield_eb_mu_comb->GetYaxis()->SetTitle("signal yield");

  hyield_eb_mu_comb->SetLabelSize(0, "X");
  hyield_eb_mu_comb->SetTitleOffset(1.5, "Y");
  hyield_eb_mu_comb->SetTitleOffset(1.3, "X");
  hyield_eb_mu_comb->Draw("e");
  hyield_eb_mu_pri->Draw("same");

  TLegend *lg = new TLegend(0.4, 0.5, 0.7, 0.64);
  lg->AddEntry(hyield_eb_mu_comb, "combine tool", "ep");
  lg->AddEntry(hyield_eb_mu_pri, "private", "ep");
  lg->SetBorderSize(0);
  lg->Draw();

  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextSize(0.04);
  tx.SetTextFont(62);
  tx.DrawLatex(0.4, 0.75, "Z#gamma #rightarrow #mu#mu#gamma, EB");

  c1->cd();
  TPad *pad12 = new TPad("pad12", "pad12", 0., 0., 1, 0.25);
  pad12->SetTopMargin(0.04);
  pad12->SetBottomMargin(0.35);
  pad12->Draw();
  pad12->cd();
  TH1F *hratio_eb_mu = (TH1F*) hyield_eb_mu_comb->Clone();
  hratio_eb_mu->Divide(hyield_eb_mu_pri);
  hratio_eb_mu->SetLineColor(4);
  hratio_eb_mu->GetYaxis()->SetTitleOffset(1.6);
  hratio_eb_mu->SetMinimum(-6);
  hratio_eb_mu->SetTitleSize(0.14, "XYZ");
  hratio_eb_mu->SetLabelSize(0.14, "XYZ");
  hratio_eb_mu->SetTitleFont(42, "XYZ");
  hratio_eb_mu->GetYaxis()->SetTitle("combine/private");
  hratio_eb_mu->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
  hratio_eb_mu->GetYaxis()->SetTitleOffset(0.55);
  hratio_eb_mu->GetXaxis()->SetTitleOffset(1.1);
  hratio_eb_mu->GetYaxis()->SetRangeUser(0.9, 1.1);
  hratio_eb_mu->GetYaxis()->SetNdivisions(8);
  hratio_eb_mu->Draw();

  TLine *l = new TLine(-1,1,1,1);
  l->SetLineWidth(2);
  l->SetLineColor(1);
  l->SetLineStyle(kDashed);
  l->Draw("same");

  TCanvas *c2 = new TCanvas("c2", "c2", 650, 650);
  c2->cd();
  TPad *pad21 = new TPad("pad21", "pad21", 0, 0.25, 1, 1);
  pad21->SetBottomMargin(0.02);
  pad21->Draw();
  pad21->cd();
  //pad21->SetLogy();
  pad21->SetFillColor(0);

  hyield_ee_mu_comb->GetYaxis()->SetTitleOffset(1.6);
  hyield_ee_mu_comb->SetTitleSize(0.05, "XYZ");
  hyield_ee_mu_comb->SetLabelSize(0.05, "XYZ");
  hyield_ee_mu_comb->SetTitleFont(42, "XYZ");
  hyield_ee_mu_comb->GetYaxis()->SetTitle("signal yield");

  hyield_ee_mu_comb->SetLabelSize(0, "X");
  hyield_ee_mu_comb->SetTitleOffset(1.5, "Y");
  hyield_ee_mu_comb->SetTitleOffset(1.3, "X");
  hyield_ee_mu_comb->Draw("e");
  hyield_ee_mu_pri->Draw("same");
  tx.DrawLatex(0.4, 0.75, "Z#gamma #rightarrow #mu#mu#gamma, EE");

  lg->Draw();

  c2->cd();
  TPad *pad22 = new TPad("pad22", "pad22", 0., 0., 1, 0.25);
  pad22->SetTopMargin(0.04);
  pad22->SetBottomMargin(0.35);
  pad22->Draw();
  pad22->cd();
  TH1F *hratio_ee_mu = (TH1F*) hyield_ee_mu_comb->Clone();
  hratio_ee_mu->Divide(hyield_ee_mu_pri);
  hratio_ee_mu->SetLineColor(4);
  hratio_ee_mu->GetYaxis()->SetTitleOffset(1.6);
  hratio_ee_mu->SetMinimum(-6);
  hratio_ee_mu->SetTitleSize(0.14, "XYZ");
  hratio_ee_mu->SetLabelSize(0.14, "XYZ");
  hratio_ee_mu->SetTitleFont(42, "XYZ");
  hratio_ee_mu->GetYaxis()->SetTitle("combine/private");
  hratio_ee_mu->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
  hratio_ee_mu->GetYaxis()->SetTitleOffset(0.55);
  hratio_ee_mu->GetXaxis()->SetTitleOffset(1.1);
  hratio_ee_mu->GetYaxis()->SetRangeUser(0.9, 1.1);
  hratio_ee_mu->GetYaxis()->SetNdivisions(8);
  hratio_ee_mu->Draw();

  c1->SaveAs("plots/fitting/ratio_signalyield_combineTool_privateCode_EB_mu.pdf");
  c2->SaveAs("plots/fitting/ratio_signalyield_combineTool_privateCode_EE_mu.pdf");


}
