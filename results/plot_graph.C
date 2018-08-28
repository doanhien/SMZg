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


void revertX(TGraphErrors *gr)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,-x,y);
   }
}

void shiftX(TGraphErrors *gr, double dx)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,x+dx,y);
   }
}

void setErrorX(TGraphErrors *gr, double dx)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double dy = gr->GetEY()[i];
      gr->SetPointError(i,dx,dy);
   }
}

void setErrorY(TGraphAsymmErrors *gr, double dy)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double dx = gr->GetEXhigh()[i];
      gr->SetPointError(i,dx,dx,dy,dy);
   }
}

void revertX(TGraphAsymmErrors *gr)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,-x,y);
   }
}

void scaleY(TGraphErrors *gr, double scale)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      double ex = gr->GetEX()[i];
      double ey = gr->GetEY()[i];
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,ex,ey*scale);
   }
}

void scaleY(TGraphAsymmErrors *gr, double scale)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      double exl = gr->GetEXlow()[i];
      double eyl = gr->GetEYlow()[i];
      double exh = gr->GetEXhigh()[i];
      double eyh = gr->GetEYhigh()[i];
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,exl,exh,eyl*scale,eyh*scale);
   }
}

TGraphAsymmErrors* theory_errors(const char* gbasename, int channel_number, int nbins, const char* pdfname)
{
   string label;
   if (string(gbasename)==string("gA1p")) label="A1p";
   else if (string(gbasename)==string("gA1m")) label="A1m";
   else if (string(gbasename)==string("gch")) label="CA";
   else if (string(gbasename)==string("gA3")) label="A3";
   else if (string(gbasename)==string("gA4")) label="A4";
   else if (string(gbasename)==string("gyieldseb")) label="yields_eb";
   else if (string(gbasename)==string("gyieldsee")) label="yields_ee";
   else if (string(gbasename)==string("gyields")) label="yields";

   string channel = (channel_number==1) ? "elec" : "muon";

   string filename = string(label) + string("_") + string(channel) + string(".dta");
   ifstream file(filename.c_str());

   double dummy;
   double *y_median = new double[nbins];
   double *A_MCFM = new double[nbins];
   double *errA_MCFM_MCFM_up = new double[nbins];
   double *errA_MCFM_MCFM_down = new double[nbins];
   double *A_MCNLO = new double[nbins];
   double *errA_MCNLO_MCFM_up = new double[nbins];
   double *errA_MCNLO_MCFM_down = new double[nbins];
   double *errA_MCNLO_MCNLO_up = new double[nbins];
   double *errA_MCNLO_MCNLO_down = new double[nbins];
   double *errA_MCNLO_scales_up = new double[nbins];
   double *errA_MCNLO_scales_down = new double[nbins];

   for (int ibin=0; ibin<nbins; ibin++)
   {
      int ibin2=ibin;
      if (label=="CA" || label=="yieldsWp" || label=="yieldsWm") ibin2 = nbins-1-ibin;
      file >> dummy >> y_median[ibin] >> A_MCFM[ibin2] >> errA_MCFM_MCFM_up[ibin2] >> errA_MCFM_MCFM_down[ibin2] >> A_MCNLO[ibin2] >> errA_MCNLO_MCFM_up[ibin2] >> errA_MCNLO_MCFM_down[ibin2] >> errA_MCNLO_MCNLO_up[ibin2] >> errA_MCNLO_MCNLO_down[ibin2] >> errA_MCNLO_scales_up[ibin2] >> errA_MCNLO_scales_down[ibin2];

      if (label == "yields_eb" || label == "yields_ee" || label == "yields")
      {
         double deta = 1.;
         // if (channel_number==1 && (ibin2==0 || ibin2==nbins-1)) deta = 0.4;
         A_MCFM[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCFM_MCFM_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCFM_MCFM_down[ibin2] *= 208.*deta*LUMI*1e-6;
         A_MCNLO[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_MCFM_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_MCFM_down[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_MCNLO_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_MCNLO_down[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_scales_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_MCNLO_scales_down[ibin2] *= 208.*deta*LUMI*1e-6;
      }
   }

   double *erry = new double[nbins];
   for (int ibin=0; ibin<nbins; ibin++) erry[ibin] = 0.25;
   if (channel_number==1) {erry[0]=0.2; erry[nbins-1]=0.2;};

   double *A = new double[nbins];
   double *dA_up = new double[nbins];
   double *dA_down = new double[nbins];

   if (string(pdfname) == "MCFM")
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_MCFM[ibin];
         dA_up[ibin] = sqrt(pow(errA_MCFM_MCFM_up[ibin],2)+pow(errA_MCNLO_scales_up[ibin],2));
         dA_down[ibin] = sqrt(pow(errA_MCFM_MCFM_down[ibin],2)+pow(errA_MCNLO_scales_down[ibin],2));
      }
   }
   else
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_MCNLO[ibin];
         dA_up[ibin] = sqrt(pow(errA_MCNLO_MCNLO_up[ibin],2)+pow(errA_MCNLO_scales_up[ibin],2)+pow(errA_MCNLO_MCFM_up[ibin],2));
         dA_down[ibin] = sqrt(pow(errA_MCNLO_MCNLO_down[ibin],2)+pow(errA_MCNLO_scales_down[ibin],2)+pow(errA_MCNLO_MCFM_down[ibin],2));
      }
   }

   for (int ibin=0; ibin<nbins; ibin++)
      cout << y_median[ibin] << " " << A[ibin] << " +" << dA_up[ibin] << " -" << dA_down[ibin] << endl;

   TGraphAsymmErrors *gr = new TGraphAsymmErrors(nbins,y_median,A,erry,erry,dA_down,dA_up);

   return gr;
}

void plot_graph(const char* fname_mcfm="graph.root", const char* fname_mcnlo="graph.root", const char *gbasename="gyields", int channel_number=1, double xmin=15, double xmax=45, double ymin=0, double ymax=1)
{

  cout << "starting plotttttttttt" << endl;
  //const char* xlabel="p_{T}^{#gamma} [GeV]";
  //const char* ylabel="d#sigma/dp_{T}^{#gamma} [fb/GeV]";

  const char* xlabel="M_{ll#gamma} [GeV]";
  const char* ylabel="d#sigma/dM_{ll#gamma} [fb/GeV]";

   gStyle->SetOptTitle(0);
   TFile *fmcfm = new TFile(fname_mcfm);
   TFile *fmcnlo = new TFile(fname_mcnlo);

   TString gname_stat = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst = Form("%s_exp_%i",gbasename,channel_number);
   cout << "get graphic" << endl;
   TGraphErrors *gexp_stat = (TGraphErrors*) fmcnlo->Get(gname_stat);
   TGraphErrors *gexp_syst = (TGraphErrors*) fmcnlo->Get(gname_syst);
   //TGraphErrors *gexp_all  = (TGraphErrors*) gname_syst->Clone();

   //cout << "number of bins" << endl;
   int nbins = gexp_stat->GetN();
   /*
   for ( int i = 0; i < nbins; i++ ) {
     float stat_err = gexp_stat->GetErrorY(i);
     float syst_err = gexp_syst->GetErrorY(i);
     gexp_all->SetPointError(i, gexp_stat->GetErrorX(i), sqrt(pow(stat_err,2)+pow(syst_err,2)));
   }
   */

   //TGraphAsymmErrors *gth_mcfm = theory_errors(gbasename, channel_number, nbins, "MCFM");
   //TGraphAsymmErrors *gth_mcnlo = theory_errors(gbasename, channel_number, nbins, "MCNLO09");
   cout << "get graphic for theoary" << endl;
   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TGraphErrors *gth_mcfm = (TGraphErrors*) fmcfm->Get(gname_th);
   TGraphErrors *gth_mcnlo = (TGraphErrors*) fmcnlo->Get(gname_th);

   //int nbins = gth_mcnlo->GetN();

   /// cout << "style for graph" << endl;
   gth_mcfm->SetLineColor(gColorMCFM);
   gth_mcfm->SetLineWidth(2);
   gth_mcfm->SetFillColor(gColorMCFM);
   gth_mcfm->SetFillStyle(3002);
   gth_mcfm->GetXaxis()->SetLimits(xmin,xmax);
   gth_mcfm->GetYaxis()->SetRangeUser(ymin,ymax);
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
   gexp_stat->GetXaxis()->SetLimits(xmin,xmax);
   gexp_stat->GetYaxis()->SetRangeUser(ymin,ymax);
   gexp_stat->GetXaxis()->SetTitle(xlabel);
   gexp_stat->GetYaxis()->SetTitle(ylabel);
   gexp_stat->GetXaxis()->SetTitleOffset(1.5);
   gexp_stat->GetYaxis()->SetTitleOffset(1.7);
   //gexp_syst->Draw("P5"); 

   cout << "get ratio between data and mc" << endl;
   TGraphErrors *gpull_mcfm = new TGraphErrors;
   TGraphErrors *gpull_nlo = new TGraphErrors;

   for (int i = 0; i < nbins; i++) {
     float pmcfm = gexp_syst->GetY()[i] / gth_mcfm->GetY()[i] ;
     float pnlo  = gexp_syst->GetY()[i] / gth_mcnlo->GetY()[i] ;
     //float err_pmcfm = gexp_syst->GetErrorY(i)/gth_mcfm->GetY()[i]; //consider theory no error
     //float err_pnlo = gexp_syst->GetErrorY(i)/gth_mcnlo->GetY()[i]; //consider MC no error

     float err_pmcfm = pmcfm*sqrt(pow(gexp_syst->GetErrorY(i)/gexp_syst->GetY()[i],2) + pow(gth_mcfm->GetErrorY(i)/gth_mcfm->GetY()[i],2)); 
     float err_pnlo = pnlo*sqrt(pow(gexp_syst->GetErrorY(i)/gth_mcnlo->GetY()[i],2) +pow(gth_mcnlo->GetErrorY(i)/gth_mcnlo->GetY()[i],2));
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
   gexp_stat->Draw("PZ");
   gexp_syst->Draw("||");

   TString plabel = "Z#gamma";
   cout << "channel: " << channel_number << " region: " << gbasename << endl;
   if ( channel_number == 1) {
     if ( string(gbasename) == "gyields" ) plabel = "Z#gamma #rightarrow ee#gamma";
     else if ( string(gbasename) == "gyieldseb" ) plabel = "Z#gamma #rightarrow ee#gamma in EB";
     else if ( string(gbasename) == "gyieldsee" ) plabel = "Z#gamma #rightarrow ee#gamma in EE";
   }
   else if (channel_number == 2) {
       if (gbasename == string("gyields") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma";
       else if (gbasename == string("gyieldseb") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma in EB";
       else if (gbasename == string("gyieldsee") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma in EE";
   }

   TH1F *dummy = new TH1F();
   //TLegend *tleg = new TLegend(0.56,0.42,0.87,0.64);
   TLegend *tleg = new TLegend(0.56,0.2,0.87,0.4);
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

   gpull_mcfm->GetYaxis()->SetTitleOffset(0.5);
   gpull_mcfm->GetYaxis()->SetLabelSize(0.14);
   gpull_mcfm->GetYaxis()->SetTitleFont(42);
   gpull_mcfm->GetYaxis()->SetTitleSize(0.15);
   gpull_mcfm->GetYaxis()->SetNdivisions(5);

   gpull_mcfm->GetXaxis()->SetLabelSize(0.14);
   gpull_mcfm->GetXaxis()->SetTitleFont(42);
   gpull_mcfm->GetXaxis()->SetTitleSize(0.15);
   gpull_mcfm->GetXaxis()->SetTitleOffset(1.1);

   gpull_mcfm->GetXaxis()->SetLimits(xmin, xmax);
   //gpull_mcfm->GetYaxis()->SetRangeUser(0.8, 1.2);
   gpull_mcfm->GetYaxis()->SetRangeUser(0.5, 1.5);
   gpull_mcfm->GetYaxis()->SetTitle("data/MC");
   gpull_mcfm->GetXaxis()->SetTitle(xlabel);
   gpull_mcfm->SetLineColor(gColorMCFM);
   gpull_mcfm->SetLineStyle(9);
   gpull_mcfm->SetLineWidth(2);
   gpull_mcfm->SetFillColor(gColorMCFM);
   gpull_mcfm->SetFillStyle(3004);
   gpull_mcfm->SetMarkerColor(gColorMCFM);
   pad2->SetGridy();
   gpull_mcfm->Draw("AP5");


   gpull_nlo->SetLineColor(gColorMCNLO);
   gpull_nlo->SetLineWidth(2);
   gpull_nlo->SetFillColor(gColorMCNLO);
   gpull_nlo->SetFillStyle(3002);
   gpull_nlo->SetMarkerColor(gColorMCNLO);
   gpull_nlo->Draw("P5");

   TLine *l = new TLine(xmin,1,xmax,1);
   l->SetLineWidth(2);
   l->SetLineColor(1);
   l->SetLineStyle(kDashed);
   l->Draw("same");


   string cat = "_ele";
   if (channel_number == 1) cat = "_ele";
   else cat = "_muon";

   //c1->SaveAs(Form("plots/%s%s_SF_mcfm_nnpdf_gencut_amcnlo_2018077_newSB_2ndFit_biasCorr.pdf", gbasename,cat.c_str()));
   //c1->SaveAs(Form("plots/%s%s_vsMllg_mcfm_nnpdf_diffSeed_gencut_amcnlo_20180628.pdf", gbasename,cat.c_str()));
   //c1->SaveAs(Form("plots/%s%s_mcfm_nnpdf_nnlo_amcnlo_newSB_2ndFit_20180629.pdf", gbasename,cat.c_str()));
   c1->SaveAs(Form("plots/%s%s_Mllg_PtFakeRate.pdf", gbasename,cat.c_str()));

}



void plot_graph_theory(const char* fname_mcfm="graph.root", const char* fname_amcnlo="graph.root", const char *gbasename="gyields", int channel_number=1, double xmin=15, double xmax=45, double ymin=0, double ymax=1)
{

  //const char* xlabel="p_{T}^{#gamma} [GeV]";
  //const char* ylabel="d#sigma/dp_{T}^{#gamma} [fb/GeV]";

  const char* xlabel="M_{ll#gamma} [GeV]";
  const char* ylabel="d#sigma/dM_{ll#gamma} [fb/GeV]";

   gStyle->SetOptTitle(0);
   TFile *fmcfm = new TFile(fname_mcfm);
   TFile *famcnlo = new TFile(fname_amcnlo);

   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TGraphErrors *gth_mcfm = (TGraphErrors*) fmcfm->Get(gname_th);
   TGraphErrors *gth_amcnlo = (TGraphErrors*) famcnlo->Get(gname_th);
   cout << "number of bin in graph" << endl;
   int nbins = gth_mcfm->GetN();

   gth_mcfm->SetLineColor(gColorMCFM);
   gth_mcfm->SetLineWidth(2);
   gth_mcfm->SetFillColor(gColorMCFM);
   gth_mcfm->SetFillStyle(3002);
   gth_mcfm->GetXaxis()->SetLimits(xmin,xmax);
   gth_mcfm->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_mcfm->GetXaxis()->SetTitle(xlabel);
   gth_mcfm->GetYaxis()->SetTitle(ylabel);


   //gth_mcfm->Print();

   gth_amcnlo->SetLineColor(gColorMCNLO);
   gth_amcnlo->SetLineStyle(9);
   gth_amcnlo->SetLineWidth(2);
   gth_amcnlo->SetFillColor(gColorMCNLO);
   gth_amcnlo->SetFillStyle(3004);

   TGraphErrors *gratio = new TGraphErrors;

   for (int i = 0; i < nbins; i++) {
     float pmcfm = gth_amcnlo->GetY()[i] / gth_mcfm->GetY()[i] ;
     float err_pmcfm = pmcfm*sqrt( pow(gth_amcnlo->GetErrorY(i)/gth_amcnlo->GetY()[i],2) + pow(gth_mcfm->GetErrorY(i)/gth_mcfm->GetY()[i],2));
     gratio->SetPoint(i, gth_mcfm->GetX()[i], pmcfm );
     gratio->SetPointError(i, gth_mcfm->GetErrorX(i), err_pmcfm );

   }

   cout << "plotting canvas" << endl;

   TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
   c1->cd();
   TPad *pad1 = new TPad("pad1", "pad1", 0., 0.25, 1.,1.);
   pad1->SetBottomMargin(0.02);
   pad1->Draw();
   pad1->cd();
   pad1->SetLogy();
   pad1->SetLogx();
   gth_mcfm->GetXaxis()->SetLabelSize(0.);
   pad1->SetLogy();
   gth_mcfm->GetXaxis()->SetLimits(xmin,xmax);
   gth_mcfm->GetXaxis()->SetTickLength(0.05);
   gth_mcfm->GetXaxis()->SetMoreLogLabels();
   gth_mcfm->GetXaxis()->SetNoExponent();
   gth_mcfm->Draw("A5");
   gth_amcnlo->Draw("P5");

   TString plabel = "Z#gamma";
   cout << "channel: " << channel_number << " region: " << gbasename << endl;
   if ( channel_number == 1) {
     if ( string(gbasename) == "gyields" ) plabel = "Z#gamma #rightarrow ee#gamma";
     else if ( string(gbasename) == "gyieldseb" ) plabel = "Z#gamma #rightarrow ee#gamma in EB";
     else if ( string(gbasename) == "gyieldsee" ) plabel = "Z#gamma #rightarrow ee#gamma in EE";
   }
   else if (channel_number == 2) {
       if (gbasename == string("gyields") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma";
       else if (gbasename == string("gyieldseb") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma in EB";
       else if (gbasename == string("gyieldsee") ) plabel = "Z#gamma #rightarrow #mu#mu#gamma in EE";
   }

   TLegend *tleg = new TLegend(0.56,0.52,0.87,0.68);
   tleg->AddEntry(gth_mcfm,"MCFM_NLO","lf");
   tleg->AddEntry(gth_amcnlo,"aMCatNLO (NNPDF)","lf");
   //tleg->AddEntry(gth_amcnlo,"MCFM_NNLO","lf");
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
   pad2->SetLogx();

   gratio->GetXaxis()->SetTickLength(0.05);
   gratio->GetXaxis()->SetMoreLogLabels();
   gratio->GetXaxis()->SetNoExponent();

   gratio->GetYaxis()->SetTitleOffset(0.5);
   gratio->GetYaxis()->SetLabelSize(0.14);
   gratio->GetYaxis()->SetTitleFont(42);
   gratio->GetYaxis()->SetTitleSize(0.15);
   gratio->GetYaxis()->SetNdivisions(5);

   gratio->GetXaxis()->SetLabelSize(0.14);
   gratio->GetXaxis()->SetTitleFont(42);
   gratio->GetXaxis()->SetTitleSize(0.15);
   gratio->GetXaxis()->SetTitleOffset(1.1);

   gratio->GetXaxis()->SetLimits(xmin, xmax);
   //gratio->GetYaxis()->SetRangeUser(0.8, 1.2);
   gratio->GetYaxis()->SetRangeUser(0.5, 1.5);
   //gratio->GetYaxis()->SetTitle("aMCatNLO/MCFM");
   gratio->GetYaxis()->SetTitle("NNLO/NLO");
   gratio->GetXaxis()->SetTitle(xlabel);
   gratio->SetLineColor(gColorMCNLO);
   gratio->SetLineStyle(9);
   gratio->SetLineWidth(2);
   gratio->SetFillColor(gColorMCNLO);
   gratio->SetFillStyle(3004);
   gratio->SetMarkerColor(gColorMCNLO);
   pad2->SetGridy();
   gratio->Draw("AP5");


   TLine *l = new TLine(xmin,1,xmax,1);
   l->SetLineWidth(2);
   l->SetLineColor(1);
   l->SetLineStyle(kDashed);
   //l->Draw("same");


   string cat = "_ele";
   if (channel_number == 1) cat = "_ele";
   else cat = "_muon";
   //gPad->SaveAs(Form("plots/%s%s_SF_mcfm.pdf", gbasename,cat.c_str()));
   c1->SaveAs(Form("plots/%s%s_mcfm_amcnlo.pdf", gbasename,cat.c_str()));

}

void plot_graph_2chan(const char* fname_mcnlo_el ="graph_mcfm_1.root", const char* fname_mcnlo_mu="graph_mcnlo_1.root", const char *gbasename="gyields", const char* xlabel="#eta_{lab}", const char* ylabel="asymetry", double xmin=15, double xmax=45, double ymin=0, double ymax=-1)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(10);
   TFile *fmcnlo_el = new TFile(fname_mcnlo_el);
   TFile *fmcnlo_mu = new TFile(fname_mcnlo_mu);

   int channel_number = 1;

   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat1 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst1 = Form("%s_exp_%i",gbasename,channel_number);

   channel_number=2;
   TString gname_th2 = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) fmcnlo_el->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) fmcnlo_el->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) fmcnlo_el->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) fmcnlo_el->Get(gname_syst2);

   int nbins = gexp_stat1->GetN();
   TGraphErrors *gth_mcnlo_el = (TGraphErrors*) fmcnlo_el->Get(gname_th);
   TGraphErrors *gth_mcnlo_mu = (TGraphErrors*) fmcnlo_el->Get(gname_th2);


   TCanvas *c1 = new TCanvas("c1", "c1", 650, 650);
   c1->cd();
   c1->SetLogy();
   //TPad *pad1 = new TPad("pad1", "pad1", 0., 0.25, 1.,1.);
   //pad1->SetBottomMargin(0.02);
   //pad1->Draw();
   //pad1->cd();
   //pad1->SetLogy();
   gth_mcnlo_el->SetLineColor(gColorMCFM);
   gth_mcnlo_el->SetLineWidth(2);
   gth_mcnlo_el->SetFillColor(gColorMCFM);
   gth_mcnlo_el->SetFillStyle(3002);
   if (xmin<xmax) gth_mcnlo_el->GetXaxis()->SetLimits(xmin,xmax);
   if (ymin<ymax) gth_mcnlo_el->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_mcnlo_el->GetXaxis()->SetTitle(xlabel);
   gth_mcnlo_el->GetYaxis()->SetTitle(ylabel);
   gth_mcnlo_el->GetXaxis()->SetLabelSize(0);
   //gth_mcnlo_el->Draw("AL3");
   gth_mcnlo_el->Draw("A5");

   gth_mcnlo_el->Print();
   cout << "\n" << endl;
   gth_mcnlo_mu->Print();

   gth_mcnlo_mu->SetLineColor(gColorMCNLO);
   gth_mcnlo_mu->SetLineStyle(9);
   gth_mcnlo_mu->SetLineWidth(2);
   gth_mcnlo_mu->SetFillColor(gColorMCNLO);
   gth_mcnlo_mu->SetFillStyle(3003);
   //gth_mcnlo_mu->Draw("PZ");

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
   //gth_mcnlo_el->Draw("LP");

   cout << "\n" << "observed:" << endl;
   gexp_stat1->Print();
   cout << "\n" << endl;
   gexp_stat2->Print(); 

   gth_mcnlo_mu->SetLineColor(kOrange-3);
   gth_mcnlo_mu->SetMarkerColor(kOrange-3);
   gth_mcnlo_mu->SetFillColor(kOrange);
   gth_mcnlo_mu->SetMarkerStyle(33);
   //gth_mcnlo_mu->SetLineWidth(2);
   gth_mcnlo_mu->SetMarkerSize(2);
   gth_mcnlo_mu->SetFillStyle(3001);
   gth_mcnlo_mu->Draw("5");

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
   

   //TLegend *tleg = new TLegend(0.44,0.5,0.76,0.8);
   TLegend *tleg = new TLegend(0.6,0.5,0.86,0.85);
   tleg->AddEntry(gexp_stat1,channeltext1,"p");
   tleg->AddEntry(gexp_stat2,channeltext2,"p");
   tleg->AddEntry(gth_mcnlo_el,"MC ee#gamma","fl");
   tleg->AddEntry(gth_mcnlo_mu,"MC #mu#mu#gamma","fl");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();

   /*
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
   */

   c1->SaveAs(Form("plots/%s_phopt_2chan_rhocorr.pdf", gbasename));

}

void plot_graph_1file(const char* fname="graph.root", const char *gbasename="gyields", int channel_number=0, double scale=1.)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(6);
   gStyle->SetHatchesLineWidth(2);
   gStyle->SetHatchesSpacing(1.7);
   TFile *f = new TFile(fname);

   string xlabel, ylabel;
   double xmin, xmax, ymin, ymax;

   // x-axis
   xlabel = string("p_{T}^{#gamma} [GeV]");
   ylabel = string("d#sigma/dp_{T}^{#gamma} [fb/GeV]");
   xmin = 15.;
   xmax = 45.;

   ymin = 0.;
   ymax = 0.9;

   // y-axis range

   /*if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm")
   {
      ymin = 0.;
      ymax = 130.;
      if (string(gbasename)=="gyieldsp") ylabel = string("d#sigma (W^{+}#rightarrowl^{+}#nu) / d#eta_{lab} [nb]");
      else ylabel = string("d#sigma (W^{-}#rightarrowl^{-}#nu) / d#eta_{lab} [nb]");
   }
   else if (string(gbasename)=="gch")
   {
      ymin = -0.4;
      ymax = 0.4;
      ylabel = string("(N^{+}-N^{-})/(N^{+}+N^{-})");
      }*/

   bool twochan = (channel_number==0 || channel_number==3);
   if (twochan) channel_number = 1;

   TString gname_stat1 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst1 = Form("%s_exp_%i",gbasename,channel_number);
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,2);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,2);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) f->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) f->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f->Get(gname_syst2);
   // TFile *fc = new TFile("combo_muacc_20140219.root");
   int nbins = gexp_stat1->GetN();
   TGraphAsymmErrors *gth_mcfm = theory_errors(gbasename, (twochan) ? 1 : channel_number, nbins, "MCFM");
   TGraphAsymmErrors *gth_mcnlo = theory_errors(gbasename, (twochan) ? 1 : channel_number, nbins, "MCNLO");

   scaleY(gexp_stat1,scale);
   scaleY(gexp_syst1,scale);
   if (twochan)
   {
      scaleY(gexp_stat2,scale);
      scaleY(gexp_syst2,scale);
      setErrorX(gexp_stat1,0);
      setErrorX(gexp_syst1,0);
      setErrorX(gexp_stat2,0);
      setErrorX(gexp_syst2,0);
      shiftX(gexp_stat1,0.1);
      shiftX(gexp_syst1,0.1);
      shiftX(gexp_stat2,0.1);
      shiftX(gexp_syst2,0.1);
   }
   else
   {
      scaleY(gth_mcfm,scale);
      scaleY(gth_mcnlo,scale);
   }

   TCanvas *c = new TCanvas("c", "c", 650, 500);
   c->cd();
   c->SetLogy();
   gth_mcfm->SetLineColor(gColorMCFM);
   gth_mcfm->SetMarkerSize(0);
   gth_mcfm->SetLineWidth(2);
   gth_mcfm->SetFillColor(gColorMCFM_fill);
   // gth_mcfm->SetFillStyle(3002);
   gth_mcfm->SetFillStyle(1001);
   // gth_mcfm->Draw("AL3");
   gth_mcfm->Draw("A2");
   gth_mcfm->GetXaxis()->SetRangeUser(xmin,xmax);
   gth_mcfm->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_mcfm->GetXaxis()->SetTitle(xlabel.c_str());
   gth_mcfm->GetYaxis()->SetTitle(ylabel.c_str());
   // gth_mcfm->Draw("AL3");
   gth_mcfm->Draw("A2");
   TGraphAsymmErrors *gth_mcfm2 = new TGraphAsymmErrors(*gth_mcfm);
   setErrorY(gth_mcfm2,0);
   gth_mcfm2->Draw("Z");

   gth_mcnlo->SetLineColor(gColorMCNLO);
   gth_mcnlo->SetLineStyle(9);
   gth_mcnlo->SetLineWidth(2);
   gth_mcnlo->SetFillColor(gColorMCNLO);
   gth_mcnlo->SetFillStyle(3345);
   // gth_mcnlo->Draw("L3");
   gth_mcnlo->Draw("2");

   gexp_syst1->SetLineColor((twochan) ? gMyColor1 : gMyColor);
   gexp_syst1->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
   gexp_syst1->SetMarkerStyle(gMyMarker1);
   gexp_syst1->SetFillStyle(0);
   gexp_syst1->SetLineWidth(2);
   gexp_syst1->Draw("||");

   gexp_stat1->SetLineColor((twochan) ? gMyColor1 : gMyColor);
   gexp_stat1->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
   gexp_stat1->SetMarkerStyle(gMyMarker1);
   gexp_stat1->SetLineWidth(2);
   gexp_stat1->GetXaxis()->SetRangeUser(xmin,xmax);
   gexp_stat1->GetYaxis()->SetRangeUser(ymin,ymax);
   gexp_stat1->GetXaxis()->SetTitle(xlabel.c_str());
   gexp_stat1->GetYaxis()->SetTitle(ylabel.c_str());
   gexp_stat1->Draw("APZ");
   gth_mcfm->Draw("PZ2");


   if (twochan)
   {
      gexp_syst2->SetLineColor((twochan) ? gMyColor2 : gMyColor);
      gexp_syst2->SetMarkerColor((twochan) ? gMyColor2 : gMyColor);
      gexp_syst2->SetMarkerStyle(gMyMarker2);
      gexp_syst2->SetFillStyle(0);
      gexp_syst2->SetLineWidth(2);
      gexp_syst2->Draw("||");

      gexp_stat2->SetLineColor((twochan) ? gMyColor2 : gMyColor);
      gexp_stat2->SetMarkerColor((twochan) ? gMyColor2 : gMyColor);
      gexp_stat2->SetMarkerStyle(gMyMarker2);
      gexp_stat2->SetLineWidth(2);
      gexp_stat2->Draw("PZ");
   }

   double xl=0.52,yl=0.4,dx=0.3,dy=0.28;
   TLegend *tleg = new TLegend(xl,yl,xl+dx,yl+dy);
   tleg->SetTextFont(42);
   tleg->AddEntry(gth_mcfm,"MCFM (NNPDF)","lf");
   tleg->AddEntry(gexp_syst1,"Z#gamma#rightarrowee#gamma","pe");
   tleg->AddEntry(gexp_syst2,"Z#gamma#rightarrow#mu#mu#gamma","pe");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();

   c->SaveAs(Form("%s_2chan_SF.pdf",gbasename));
}

void plot_graph_1file_all(const char* file, int chan, bool reverse_eta=true, double scale=1./34.62)
{
   plot_graph_1file(file,"gA3",chan);
   plot_graph_1file(file,"gch",chan);
   plot_graph_1file(file,"gA1p",chan);
   plot_graph_1file(file,"gA1m",chan);
   plot_graph_1file(file,"gA4",chan);
   plot_graph_1file(file,"gyieldsp",chan,scale);
   plot_graph_1file(file,"gyieldsm",chan,scale);
}
