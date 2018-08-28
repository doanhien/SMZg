#include "TFile.h"
#include "TGraph.h"

void ratio_emu() {


  gStyle->SetOptStat(0);

  //TFile *fele = new TFile("graph_ele_amcnlo_20180628.root", "read");
  //TFile *fmu = new TFile("graph_mu_amcnlo_20180628.root", "read");
  TFile *fele = new TFile("graph_Mllg_ele_amcnlo_20180628.root", "read");
  TFile *fmu = new TFile("graph_Mllg_mu_amcnlo_20180628.root", "read");

  TGraphErrors *gr_data_ele = (TGraphErrors*) fele->Get("gyields_exp_1");
  TGraphErrors *gr_mcfm_ele = (TGraphErrors*) fele->Get("gyields_th_1");

  TGraphErrors *gr_data_mu = (TGraphErrors*) fmu->Get("gyields_exp_2");
  TGraphErrors *gr_mcfm_mu = (TGraphErrors*) fmu->Get("gyields_th_2");

  TGraphErrors *gr_ratio_da = new TGraphErrors();
  TGraphErrors *gr_ratio_mcfm = new TGraphErrors();

 int  nbins = gr_data_ele->GetN();

  for (int i = 0; i < nbins; i++) {
    float pmcfm = gr_mcfm_ele->GetY()[i] / gr_mcfm_mu->GetY()[i] ;
    float pda  = gr_data_ele->GetY()[i] / gr_data_mu->GetY()[i] ;

    gr_ratio_mcfm->SetPoint(i, gr_data_ele->GetX()[i], pmcfm );
    gr_ratio_da->SetPoint(i, gr_data_ele->GetX()[i], pda );

    gr_ratio_mcfm->SetPointError(i, gr_data_ele->GetErrorX(i), 0. );
    gr_ratio_da->SetPointError(i, gr_data_ele->GetErrorX(i), 0. );
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 650, 500);
  c1->cd();
  c1->SetGridy();
  c1->SetLogx();
  gr_ratio_mcfm->GetYaxis()->SetTitle("ele/mu");
  gr_ratio_mcfm->GetXaxis()->SetTitle("p_{T}^{#gamma} (GeV)");
  //gr_ratio_mcfm->GetXaxis()->SetLimits(15, 1000);
  gr_ratio_mcfm->GetXaxis()->SetLimits(50, 1000);
  gr_ratio_mcfm->GetYaxis()->SetRangeUser(0.8, 1.2);
  gr_ratio_mcfm->Draw("APZ");
  gr_ratio_mcfm->SetMarkerColor(kBlue);
  gr_ratio_mcfm->SetLineColor(kBlue);
  //gr_ratio_mcfm->Draw("APZ");
  TLegend *leg = new TLegend(0.4, 0.7, 0.6, 0.82);
  leg->AddEntry(gr_ratio_da, "data", "p");
  leg->AddEntry(gr_ratio_mcfm, "mcfm", "p");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->Draw();

  //c1->SaveAs("plots/ratio_ele_mu_pt25.pdf");

  const int nbin_mass = 11;
  //float Mllg[nbin_mass+1] = {50, 85, 95, 110, 135, 170, 210, 270, 350, 470, 640, 1000};
  float Mllg[nbin_mass] = {67.5, 90, 102.5, 122.5, 152.5, 190, 240, 310, 410, 555, 820};
  float err_Mllg[nbin_mass] = {17.5, 5, 7.5, 12.5, 17.5, 20, 30, 40, 60, 85, 180};

  float yield_fit[nbin_mass] = {846.036, 11492.8, 2856.91, 5616.59, 4440.17, 2138.52, 1458.07, 874.396, 460.07, 155.538, 93.2661};
  float err_yield_fit[nbin_mass] = {33.4413, 112.809, 79.7147, 119.193, 95.3402, 68.2725, 55.65, 38.5097, 24.7503, 14.0809, 9.62201};
  float yield_fkrate[nbin_mass] = {404.775, 4572.61, 2825.51, 9091.24, 6378.52, 3255.11, 2053.75, 1062.8, 491.444, 166.705, 73.4726};
  float err_yield_fkrate[nbin_mass] = {3.90103, 40.3928, 28.0092, 78.5423, 53.3737, 29.8887, 21.7654, 13.3208, 7.89468, 3.67974, 2.0634};

  TFile *ffit = new TFile("graph_Mllg_ele_mcfm_nnlo_20180817.root", "read");
  TFile *ffkrate = new TFile("graph_Mllg_ele_PtFakeRate_20180819.root", "read");

  //TGraphErrors *gr_fit = (TGraphErrors*) ffit->Get("gyields_exp_1");
  //TGraphErrors *gr_fkrate = (TGraphErrors*) ffkrate->Get("gyields_exp_1");

  TGraphErrors *gr_fit = new TGraphErrors(nbin_mass, Mllg, yield_fit, err_Mllg, err_yield_fit);
  TGraphErrors *gr_fkrate = new TGraphErrors(nbin_mass, Mllg, yield_fkrate, err_Mllg, err_yield_fkrate);

  gr_fkrate->SetMarkerColor(2);
  gr_fkrate->SetMarkerStyle(21);
  gr_fkrate->SetLineColor(2);

  gr_fit->SetMarkerColor(1);
  gr_fit->SetMarkerStyle(20);
  gr_fit->SetLineColor(1);

  //get ratio 
  TGraphErrors *gr_ratio_yield = new TGraphErrors();

  //int  nbin_mass = gr_fit->GetN();

  for (int i = 0; i < nbin_mass; i++) {
    float ratio_yield = gr_fit->GetY()[i] / gr_fkrate->GetY()[i] ;
    float err_fit = gr_fit->GetErrorY(i);
    float err_fkrate = gr_fkrate->GetErrorY(i);
    float err_ratio = ratio_yield*TMath::Sqrt(TMath::Power(err_fit/gr_fit->GetY()[i],2) + TMath::Power(err_fkrate/gr_fkrate->GetY()[i],2));

    gr_ratio_yield->SetPoint(i, gr_fit->GetX()[i], ratio_yield );
    gr_ratio_yield->SetPointError(i, gr_fit->GetErrorX(i), err_ratio );

  }


  c2 = new TCanvas("c2", "c2", 650, 650);
  c2->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetLogx();
  pad1->SetFillColor(0);

  gr_fit->GetYaxis()->SetTitleOffset(1.6);
  gr_fit->GetYaxis()->SetTitleSize(0.05);
  gr_fit->GetXaxis()->SetTitleSize(0.05);
  gr_fit->GetYaxis()->SetLabelSize(0.05);
  gr_fit->GetXaxis()->SetLabelSize(0.05);
  gr_fit->GetYaxis()->SetTitleFont(42);
  gr_fit->GetXaxis()->SetTitleFont(42);

  gr_fit->GetXaxis()->SetLabelSize(0);
  gr_fit->GetYaxis()->SetTitleOffset(1.5);
  gr_fit->GetXaxis()->SetTitleOffset(1.3);

  gr_fit->GetYaxis()->SetTitle("signal yield");
  gr_fit->GetYaxis()->SetRangeUser(70, 13000);
  gr_fit->GetXaxis()->SetLimits(50, 1000);

  gr_fit->Draw("APZ");
  gr_fkrate->Draw("PZ");

  TLegend *lg = new TLegend(0.3, 0.3, 0.6, 0.42);
  lg->AddEntry(gr_fit, "from Fit", "p");
  lg->AddEntry(gr_fkrate, "from Pt Fake Rate", "p");
  lg->SetBorderSize(0);
  lg->Draw();

  c2->cd();
  TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();
  pad2->SetGridy();
  pad2->SetLogx();
  gr_ratio_yield->GetYaxis()->SetTitleOffset(1.6);
  gr_ratio_yield->SetMinimum(-6);
  gr_ratio_yield->GetYaxis()->SetTitleSize(0.14);
  gr_ratio_yield->GetYaxis()->SetLabelSize(0.14);
  gr_ratio_yield->GetYaxis()->SetTitleFont(42);
  gr_ratio_yield->GetXaxis()->SetTitleSize(0.14);
  gr_ratio_yield->GetXaxis()->SetLabelSize(0.14);
  gr_ratio_yield->GetXaxis()->SetTitleFont(42);

  gr_ratio_yield->GetYaxis()->SetTitle("fit/FakeRate");
  gr_ratio_yield->GetXaxis()->SetTitle("M_{ll#gamma} [GeV]");
  gr_ratio_yield->GetYaxis()->SetTitleOffset(0.55);
  gr_ratio_yield->GetXaxis()->SetTitleOffset(1.1);
  gr_ratio_yield->GetYaxis()->SetRangeUser(0., 2.0);
  gr_ratio_yield->GetYaxis()->SetLimits(50, 1000);
  gr_ratio_yield->GetYaxis()->SetNdivisions(8);
  gr_ratio_yield->GetXaxis()->SetMoreLogLabels();

  gr_ratio_yield->Draw("APZ");


  c2->SaveAs("plots/ratio_Mllg_FromFit_FromRatevsPt.pdf");




}



