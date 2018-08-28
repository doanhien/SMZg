#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include <fstream>
#include <string>
#include "lumierror.h"
#include "CMS_lumi.C"

// EColor gMyColor = kViolet;
EColor gMyColor = kRed;
EColor gMyColor1 = kRed;
EColor gMyColor2 = kBlue;
EColor gColorMCFM = kBlue;
EColor gColorMCFM_fill = kBlue;
EColor gColorMCNLO = kGreen;
int gMyMarker = 22;
int gMyMarker1 = 20;
int gMyMarker2 = 21;

TString channeltext = "Z#gamma #rightarrow ee#gamma";
const char *channeltext_p = "Z#gamma #rightarrow ee#gamma in EB";
const char *channeltext_m = "Z#gamma #rightarrow ee#gamma in EE";
const char *channeltext1 = "Data Z#gamma #rightarrow ee#gamma";
const char *channeltext1_p = "Z#gamma #rightarrow ee#gamma in EB";
const char *channeltext1_m = "Z#gamma #rightarrow ee#gamma in EE";
const char *channeltext2 = "Data Z#gamma #rightarrow #mu#mu#gamma";
const char *channeltext2_p = "Z#gamma #rightarrow #mu#mu#gamma in EB";
const char *channeltext2_m = "Z#gamma #rightarrow #mu#mu#gamma in EE";


void plot_graph(const char* fname_mcfm="graph.root", const char* fname_mcnlo="graph.root", const char *gbasename="gyields", const char *cat="pt", int channel_number=1, double xmin=15, double xmax=45, double ymin=0, double ymax=1)
{

  cout << "starting plotttttttttt" << endl;

  string xlabel, ylabel;
  if ( string(cat)=="pt" ) {
    xlabel="p_{T}^{#gamma} [GeV]";
    ylabel="d#sigma/dp_{T}^{#gamma} [fb/GeV]";}
  else {
    xlabel="M_{ll#gamma} [GeV]";
    ylabel="d#sigma/dM_{ll#gamma} [fb/GeV]";
  }


   gStyle->SetOptTitle(0);
   TFile *fmcfm = new TFile(fname_mcfm);
   TFile *fmcnlo = new TFile(fname_mcnlo);

   TString gname_stat = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gexp_stat = (TGraphErrors*) fmcnlo->Get(gname_stat);
   TGraphErrors *gexp_syst = (TGraphErrors*) fmcnlo->Get(gname_syst);

   int nbins = gexp_stat->GetN();
   //TGraphAsymmErrors *gth_mcfm = theory_errors(gbasename, channel_number, nbins, "MCFM");
   //TGraphAsymmErrors *gth_mcnlo = theory_errors(gbasename, channel_number, nbins, "MCNLO09");

   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TGraphErrors *gth_mcfm = (TGraphErrors*) fmcfm->Get(gname_th);
   TGraphErrors *gth_mcnlo = (TGraphErrors*) fmcnlo->Get(gname_th);


   //convert from tgraph to histogram
   float xbin[nbins+1];;
   if ( string(cat) == "pt") {
     xbin[0] = 15, xbin[1] =  20, xbin[2] = 25, xbin[3] = 30, xbin[4] = 35, xbin[5] = 45;
     xbin[6] = 55, xbin[7] = 65, xbin[8] = 75, xbin[9] = 85, xbin[10] = 95;
     xbin[11] = 120, xbin[12] = 1000;
   }
   else {
     xbin[0] = 50, xbin[1] = 85, xbin[2] = 95 ,xbin[3] = 110, xbin[4] = 135;
     xbin[5] = 170, xbin[6] = 210, xbin[7] = 270, xbin[8] = 350, xbin[9] = 470;
     xbin[10] = 640, xbin[11] = 1000;
   }

   TH1F *hexp_stat = new TH1F("hexp_stat", "hexp_stat", nbins, xbin);
   TH1F *hexp_syst = new TH1F("hexp_syst", "hexp_syst", nbins, xbin);
   TH1F *hth_mcfm = new TH1F("hth_mcfm", "hth_mcfm", nbins, xbin);
   TH1F *hth_mcnlo = new TH1F("hth_mcnlo", "hth_mcnlo", nbins, xbin);

   for (int i = 0; i < nbins; i++) {
     double xstat,ystat;
     double xsyst,ysyst;
     double xmcfm,ymcfm;
     double xmcnlo,ymcnlo;

     gexp_stat->GetPoint(i, x, y); 
     hexp_stat->Fill(x,y)



   //int nbins = 12;

   /// cout << "style for graph" << endl;
   gth_mcfm->SetLineColor(gColorMCFM);
   gth_mcfm->SetLineWidth(2);
   gth_mcfm->SetFillColor(gColorMCFM);
   gth_mcfm->SetFillStyle(3002);
   if (xmin<xmax) gth_mcfm->GetXaxis()->SetLimits(xmin,xmax);
   if (ymin<ymax) gth_mcfm->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_mcfm->GetXaxis()->SetTitle(xlabel);
   gth_mcfm->GetYaxis()->SetTitle(ylabel);


   //gth_mcfm->Print();

   gth_mcnlo->SetLineColor(gColorMCNLO);
   gth_mcnlo->SetLineStyle(9);
   gth_mcnlo->SetLineWidth(2);
   gth_mcnlo->SetFillColor(gColorMCNLO);
   gth_mcnlo->SetFillStyle(3004);

   gexp_syst->SetLineColor(gMyColor);
   gexp_syst->SetMarkerColor(gMyColor);
   gexp_syst->SetMarkerStyle(8);
   gexp_syst->SetFillStyle(0);
   gexp_syst->SetLineWidth(2);
   //gexp_syst->Draw("P5");

   gexp_stat->SetLineColor(gMyColor);
   gexp_stat->SetMarkerColor(gMyColor);
   gexp_stat->SetMarkerStyle(8);
   gexp_stat->SetLineWidth(2);
   if (xmin<xmax) gexp_stat->GetXaxis()->SetLimits(xmin,xmax);
   if (ymin<ymax) gexp_stat->GetYaxis()->SetRangeUser(ymin,ymax);
   gexp_stat->GetXaxis()->SetTitle(xlabel);
   gexp_stat->GetYaxis()->SetTitle(ylabel);
   gexp_stat->GetXaxis()->SetTitleOffset(1.5);
   gexp_stat->GetYaxis()->SetTitleOffset(1.7);
   //gexp_syst->Draw("P5"); 

   cout << "get ratio between data and mc" << endl;
   TGraphErrors *gpull_mcfm = new TGraphErrors;
   TGraphErrors *gpull_nlo = new TGraphErrors;

   for (int i = 0; i < nbins; i++) {
     //float pmcfm = gexp_stat->GetY()[i] / gth_mcfm->GetY()[i] ;
     //float pnlo  = gexp_stat->GetY()[i] / gth_mcnlo->GetY()[i] ;
     ////float err_pmcfm = gexp_stat->GetErrorY(i)/gth_mcfm->GetY()[i]; 
     //float err_pmcfm = pmcfm*sqrt(pow(gexp_stat->GetErrorY(i)/gexp_stat->GetY()[i],2) + pow(gth_mcfm->GetErrorY(i)/gth_mcfm->GetY()[i],2));
     //float err_pnlo = gexp_stat->GetErrorY(i)/gth_mcnlo->GetY()[i]; //consider MC no error

     float pmcfm = gexp_syst->GetY()[i] / gth_mcfm->GetY()[i] ;
     float pnlo  = gexp_syst->GetY()[i] / gth_mcnlo->GetY()[i] ;
     //float err_pmcfm = gexp_syst->GetErrorY(i)/gth_mcfm->GetY()[i]; 
     float err_pmcfm = pmcfm*sqrt(pow(gexp_syst->GetErrorY(i)/gexp_syst->GetY()[i],2) + pow(gth_mcfm->GetErrorY(i)/gth_mcfm->GetY()[i],2)); 
     float err_pnlo = gexp_syst->GetErrorY(i)/gth_mcnlo->GetY()[i]; //consider MC no error
     gpull_mcfm->SetPoint(i, gth_mcfm->GetX()[i], pmcfm );
     gpull_mcfm->SetPointError(i, gth_mcfm->GetErrorX(i), err_pmcfm );

     gpull_nlo->SetPoint(i, gth_mcfm->GetX()[i], pnlo );
     gpull_nlo->SetPointError(i, gth_mcfm->GetErrorX(i), err_pnlo );

     }

   cout << "plotting" << endl;
   TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0., 0.25, 1.,1.);
   pad1->SetBottomMargin(0.02);
   pad1->Draw();
   pad1->cd();
   gth_mcfm->GetXaxis()->SetLabelSize(0.);
   pad1->SetLogy();
   gth_mcfm->Draw("A5");
   gth_mcnlo->Draw("P5");
   //gexp_stat->Draw("PZ");
   gexp_syst->Draw("PZ");

   TString plabel = "Z#gamma";
   cout << "channel: " << channel_number << " region: " << gbasename << endl;
   if ( channel_number == 1) {
     if ( string(gbasename) == "gyields" ) plabel = "Z#gamma #rightarrow ee#gamma";
   }
   else if (channel_number == 2) {
       if (gbasename == string("gyields") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma";
   }

   TH1F *dummy = new TH1F();
   //TLegend *tleg = new TLegend(0.56,0.42,0.87,0.64);
   TLegend *tleg = new TLegend(0.56,0.2,0.87,0.4);
   //tleg->AddEntry(dummy,"CMS Preliminary","");
   //tleg->AddEntry(dummy,"pp 13 TeV (2016)","");
   //tleg->AddEntry(dummy,"L = 35.9 fb^{-1}","");
   tleg->AddEntry(gexp_stat,plabel,"p");
   //tleg->AddEntry(gth_mcfm,"MCFM_CT14","lf");
   tleg->AddEntry(gth_mcfm,"MCFM_NPDF","lf");
   tleg->AddEntry(gth_mcnlo,"aMCatNLO_NNPDF","lf");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->SetTextSize(0.05);
   tleg->Draw();

   CMS_lumi(pad1, 4, 0);

   c1->cd();
   TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
   pad2->SetTopMargin(0.04);
   pad2->SetBottomMargin(0.35);
   pad2->Draw();
   pad2->cd();

   gpull_nlo->GetYaxis()->SetTitleOffset(0.5);
   gpull_nlo->GetYaxis()->SetLabelSize(0.14);
   gpull_nlo->GetYaxis()->SetTitleFont(42);
   gpull_nlo->GetYaxis()->SetTitleSize(0.15);
   gpull_nlo->GetYaxis()->SetNdivisions(5);

   gpull_nlo->GetXaxis()->SetLabelSize(0.14);
   gpull_nlo->GetXaxis()->SetTitleFont(42);
   gpull_nlo->GetXaxis()->SetTitleSize(0.15);
   gpull_nlo->GetXaxis()->SetTitleOffset(1.1);

   gpull_nlo->GetXaxis()->SetLimits(xmin, xmax);
   //gpull_nlo->GetYaxis()->SetRangeUser(0.8, 1.2);
   gpull_nlo->GetYaxis()->SetRangeUser(0.5, 1.5);
   gpull_nlo->GetYaxis()->SetTitle("data/MC");
   gpull_nlo->GetXaxis()->SetTitle(xlabel);
   gpull_nlo->SetLineColor(gColorMCNLO);
   gpull_nlo->SetLineStyle(9);
   gpull_nlo->SetLineWidth(2);
   gpull_nlo->SetFillColor(gColorMCNLO);
   gpull_nlo->SetFillStyle(3004);
   gpull_nlo->SetMarkerColor(gColorMCNLO);
   pad2->SetGridy();
   gpull_nlo->Draw("AP5");


   gpull_mcfm->SetLineColor(gColorMCFM);
   gpull_mcfm->SetLineWidth(2);
   gpull_mcfm->SetFillColor(gColorMCFM);
   gpull_mcfm->SetFillStyle(3002);
   gpull_mcfm->SetMarkerColor(gColorMCFM);
   gpull_mcfm->Draw("PZ");

   TLine *l = new TLine(xmin,1,xmax,1);
   l->SetLineWidth(2);
   l->SetLineColor(1);
   l->SetLineStyle(kDashed);
   l->Draw("same");


   string cat = "_ele";
   if (channel_number == 1) cat = "_ele";
   else cat = "_muon";

   //c1->SaveAs(Form("plots/%s%s_SF_mcfm_nnpdf_diffSeed_gencut_amcnlo_20180628.pdf", gbasename,cat.c_str()));
   c1->SaveAs(Form("plots/%s%s_vsMllg_mcfm_nnpdf_diffSeed_gencut_amcnlo_20180628.pdf", gbasename,cat.c_str()));

}


void plot_graph_2chan(const char* fname_mcfm1="graph_mcfm_1.root", const char* fname_eps1="graph_eps_1.root", const char* fname_2="graph_2.root", const char *gbasename="gyields", const char* xlabel="#eta_{lab}", const char* ylabel="asymetry", double xmin=15, double xmax=45, double ymin=0, double ymax=-1)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(10);
   TFile *fmcfm1 = new TFile(fname_mcfm1);
   TFile *feps1 = new TFile(fname_eps1);
   TFile *f2 = new TFile(fname_2);

   int channel_number = 1;


   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat1 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst1 = Form("%s_exp_%i",gbasename,channel_number);
   channel_number=2;
   TString gname_th2 = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) fmcfm1->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) fmcfm1->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f2->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f2->Get(gname_syst2);
   int nbins = gexp_stat1->GetN();
   TGraphErrors *gth_mcfm = (TGraphErrors*) fmcfm1->Get(gname_th);
   TGraphErrors *gth_mcfm2 = (TGraphErrors*) fmcfm1->Get(gname_th2);
   TGraphErrors *gth_eps = (TGraphErrors*) feps1->Get(gname_th);

   
   TGraphErrors *gpull_ele = new TGraphErrors;
   TGraphErrors *gpull_mu = new TGraphErrors;

   //cout << "number of x points: " << gexp_stat2->GetN() << endl;
   for (int i = 0; i < nbins; i++) {
     //float pele = (gexp_stat1->GetY()[i] - gth_mcfm->GetY()[i]) / gexp_stat1->GetErrorY(i);
     float pele = gexp_stat1->GetY()[i] / gth_mcfm->GetY()[i] ;
     float pmu  = gexp_stat2->GetY()[i] / gth_mcfm2->GetY()[i] ;
     //float err_pele = pow(gexp_stat1->GetErrorY(i)/gexp_stat1->GetY()[i],2) + pow(gth_mcfm->GetErrorY(i)/gth_mcfm->GetY()[i],2);
     //err_pele = sqrt(err_pele);
     //err_pele *= pele;
     float err_pele = gexp_stat1->GetErrorY(i)/gth_mcfm->GetY()[i]; //consider MC no error
     float err_pmu = gexp_stat2->GetErrorY(i)/gth_mcfm2->GetY()[i];
     //cout << pele << "+/-" << err_pele << endl;
     gpull_ele->SetPoint(i, gth_mcfm->GetX()[i], pele );
     gpull_ele->SetPointError(i, gth_mcfm->GetErrorX(i), err_pele );

     gpull_mu->SetPoint(i, gth_mcfm->GetX()[i], pmu );
     gpull_mu->SetPointError(i, gth_mcfm->GetErrorX(i), err_pmu );

     }


   TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0., 0.25, 1.,1.);
   pad1->SetBottomMargin(0.02);
   pad1->Draw();
   pad1->cd();
   pad1->SetLogy();
   //gth_mcfm->SetLineColor(gColorMCFM);
   //gth_mcfm->SetLineWidth(2);
   //gth_mcfm->SetFillColor(gColorMCFM);
   gth_mcfm->SetFillStyle(3002);
   //gth_mcfm->SetMarkerStyle(24);
   //gth_mcfm->SetMarkerSize(1);
   gth_mcfm->SetLineColor(kGreen+1);
   gth_mcfm->SetFillColor(kGreen+1);
   gth_mcfm->SetMarkerColor(kGreen+1);
   if (xmin<xmax) gth_mcfm->GetXaxis()->SetLimits(xmin,xmax);
   if (ymin<ymax) gth_mcfm->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_mcfm->GetXaxis()->SetTitle(xlabel);
   gth_mcfm->GetYaxis()->SetTitle(ylabel);
   gth_mcfm->GetXaxis()->SetLabelSize(0);
   //gth_mcfm->Draw("AL3");
   gth_mcfm->Draw("A5");

   gth_mcfm->Print();
   cout << "\n" << endl;
   gth_mcfm2->Print();

   gth_eps->SetLineColor(gColorMCNLO);
   gth_eps->SetLineStyle(9);
   gth_eps->SetLineWidth(2);
   gth_eps->SetFillColor(gColorMCNLO);
   gth_eps->SetFillStyle(3003);
   //gth_eps->Draw("PZ");

   gexp_syst1->SetLineColor(gMyColor1);
   gexp_syst1->SetMarkerColor(gMyColor1);
   gexp_syst1->SetMarkerStyle(8);
   gexp_syst1->SetFillStyle(0);
   gexp_syst1->SetLineWidth(2);
   //gexp_syst1->Draw("P5");

   //gPad->SetLogy();
   gexp_stat1->SetLineColor(gMyColor1);
   gexp_stat1->SetMarkerColor(gMyColor1);
   gexp_stat1->SetMarkerStyle(8);
   gexp_stat1->SetLineWidth(2);
   gexp_stat1->GetXaxis()->SetLimits(xmin,xmax);
   gexp_stat1->GetYaxis()->SetRangeUser(ymin,ymax);
   gexp_stat1->GetXaxis()->SetTitle(xlabel);
   gexp_stat1->GetYaxis()->SetTitle(ylabel);
   gexp_stat1->GetYaxis()->SetTitleOffset(1.6);
   gexp_stat1->GetXaxis()->SetTitleOffset(1.4);
   gexp_stat1->Draw("PZ");
   //gth_mcfm->Draw("LP");

   cout << "\n" << "observed:" << endl;
   gexp_stat1->Print();
   cout << "\n" << endl;
   gexp_stat2->Print(); 

   gth_mcfm2->SetLineColor(kOrange-3);
   gth_mcfm2->SetMarkerColor(kOrange-3);
   gth_mcfm2->SetFillColor(kOrange);
   gth_mcfm2->SetMarkerStyle(33);
   //gth_mcfm2->SetLineWidth(2);
   gth_mcfm2->SetMarkerSize(2);
   gth_mcfm2->SetFillStyle(3001);
   gth_mcfm2->Draw("5");

   //gexp_syst1->Print();
   //gexp_stat1->Print();

   gexp_syst2->SetLineColor(gMyColor2);
   gexp_syst2->SetMarkerColor(gMyColor2);
   gexp_syst2->SetMarkerStyle(8);
   gexp_syst2->SetFillStyle(0);
   gexp_syst2->SetLineWidth(2);
   //gexp_syst2->Draw("P5");

   gexp_stat2->SetLineColor(gMyColor2);
   gexp_stat2->SetMarkerColor(gMyColor2);
   gexp_stat2->SetMarkerStyle(22);
   gexp_stat2->SetLineWidth(2);
   gexp_stat2->Draw("PZ");
   //gexp_stat2->Print();
   

   TH1F *dummy = new TH1F();
   //TLegend *tleg = new TLegend(0.44,0.5,0.76,0.8);
   TLegend *tleg = new TLegend(0.6,0.5,0.86,0.85);
   tleg->AddEntry(dummy,"CMS Preliminary","");
   tleg->AddEntry(dummy,"pp #sqrt{s_{NN}} = 13 TeV","");
   tleg->AddEntry(dummy,"L = 35.9 fb^{-1}","");
   tleg->AddEntry(gexp_stat1,channeltext1,"p");
   tleg->AddEntry(gexp_stat2,channeltext2,"p");
   tleg->AddEntry(gth_mcfm,"MC ee#gamma","fl");
   tleg->AddEntry(gth_mcfm2,"MC #mu#mu#gamma","fl");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();


   c1->cd();
   TPad *pad2 = new TPad("pad2","pad2",0, 0, 1, 0.25);
   pad2->SetTopMargin(0.04);
   pad2->SetBottomMargin(0.35);
   pad2->Draw();
   pad2->cd();

   gpull_ele->GetYaxis()->SetTitleOffset(0.5);
   gpull_ele->GetYaxis()->SetLabelSize(0.14);
   gpull_ele->GetYaxis()->SetTitleFont(42);
   gpull_ele->GetYaxis()->SetTitleSize(0.15);
   gpull_ele->GetYaxis()->SetNdivisions(5);

   gpull_ele->GetXaxis()->SetLabelSize(0.14);
   gpull_ele->GetXaxis()->SetTitleFont(42);
   gpull_ele->GetXaxis()->SetTitleSize(0.15);
   gpull_ele->GetXaxis()->SetTitleOffset(1.1);

   gpull_ele->GetXaxis()->SetLimits(xmin, xmax);
   gpull_ele->GetYaxis()->SetRangeUser(0.8, 1.2);
   //gpull_ele->GetYaxis()->SetRangeUser(0.5, 1.5);
   gpull_ele->GetYaxis()->SetTitle("data/MC");
   gpull_ele->GetXaxis()->SetTitle(xlabel);
   gpull_ele->SetMarkerStyle(20);
   gpull_ele->SetMarkerColor(2);
   gpull_ele->SetLineColor(2);
   gpull_ele->Draw("APZ");

   gpull_mu->SetMarkerStyle(22);
   gpull_mu->SetMarkerColor(4);
   gpull_mu->SetLineColor(4);
   gpull_mu->Draw("PZ");

   TLine *l = new TLine(xmin,1,xmax,1);
   l->SetLineWidth(2);
   l->SetLineColor(1);
   l->SetLineStyle(kDashed);
   l->Draw("same");

   TLegend *le = new TLegend(0.35, 0.75, 0.5, 0.82);
   le->AddEntry(gpull_ele, "Z#gamma #rightarrow ee#gamma", "p");
   le->SetBorderSize(0);
   le->SetTextSize(0.1);
   //le->Draw();

   TLegend *lm = new TLegend(0.5, 0.75, 0.75, 0.82);
   lm->AddEntry(gpull_mu, "Z#gamma #rightarrow #mu#mu#gamma", "p");
   lm->SetBorderSize(0);
   lm->SetTextSize(0.1);
   //lm->Draw();


   c1->SaveAs(Form("plots/%s_phopt_2chan_rhocorr.pdf", gbasename));

}

