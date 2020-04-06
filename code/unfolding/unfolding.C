#include <iostream>
#include <cmath>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TROOT.h>
#include <TText.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH1.h>
#include <TF1.h>
#include <TFitter.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>
#include <TMatrixDSymEigen.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include "TUnfoldDensity.h"
#include "TUnfoldSys.h"

using namespace std;

TH2 *gHistInvEMatrix;

void unfolding( bool data = true, bool ele = true, bool barrel = true)
{

  TString finame = "inputUnfold/SB_EE7to13_EE6to14/main/ZgToLLg_5f_responseMatrix_";
  //TString suf = "isPVGood_gencut.root";
  TString suf = "isPVGood_ele.root";

  finame += suf;

  TFile* fresp = TFile::Open(finame);

  TFile *finput = new TFile("inputUnfold/SB_EE7to13_EE6to14/main/inputUnfold_isPVGood_Madgraph_data.root", "read");

  cout << "input file: " << finame.Data() << endl;

  TH1::SetDefaultSumw2();
  gStyle->SetOptFit(1111);

  Double_t const lumiData=35900;
  Double_t const lumiMC  =35900;


  TString matrix_name;
  if (ele && barrel) matrix_name = "hrec_gen_ele_EB";
  else if (ele && !barrel) matrix_name = "hrec_gen_ele_EE";
  else if (!ele && barrel) matrix_name = "hrec_gen_mu_EB";
  else if (!ele && !barrel) matrix_name = "hrec_gen_mu_EE";

  TString matrix_name_mllg;
  if (ele && barrel) matrix_name_mllg = "hrec_gen_mllg_EB_ele";
  else if (ele && !barrel) matrix_name_mllg = "hrec_gen_mllg_EE_ele";
  else if (!ele && barrel) matrix_name_mllg = "hrec_gen_mllg_EB_mu";
  else matrix_name_mllg = "hrec_gen_mllg_EE_mu";

  TString matrix_name_ptllg;
  if (ele && barrel) matrix_name_ptllg = "hrec_gen_ptllg_EB_ele";
  else if (ele && !barrel) matrix_name_ptllg = "hrec_gen_ptllg_EE_ele";
  else if (!ele && barrel) matrix_name_ptllg = "hrec_gen_ptllg_EB_mu";
  else matrix_name_ptllg = "hrec_gen_ptllg_EE_mu";

  TString matrix_name_njet;
  if (ele) matrix_name_njet = "hrec_gen_ele_njet";
  else matrix_name_njet = "hrec_gen_mu_njet";

  TH2D *histUnfoldMatrix = (TH2D*) fresp->Get(matrix_name);
  TH2D *histUnfoldMatrix_mllg = (TH2D*) fresp->Get(matrix_name_mllg);
  TH2D *histUnfoldMatrix_ptllg = (TH2D*) fresp->Get(matrix_name_ptllg);
  TH2D *histUnfoldMatrix_njet = (TH2D*) fresp->Get(matrix_name_njet);

  cout << "done reading matrix" << endl;

  TH1D *histUnfoldInput, *histUnfoldInput_mllg, *histUnfoldInput_njet;
  TH1D *histUnfoldInput_ptllg;
  TH1D *hgen, *hgen_mllg, *hgen_njet;
  TH1D *hgen_ptllg;

  TString name_InputHist, name_InputHist_mllg, name_InputHist_njet;
  TString name_InputHist_ptllg;
  TString name_hgen, name_hgen_mllg, name_hgen_njet;
  TString name_hgen_ptllg;

  if (data) {
    if (ele && barrel) {
      name_InputHist = "hda_ele_EB"; name_InputHist_mllg = "hda_mllg_EB_ele";
      name_InputHist_ptllg = "hda_ptllg_EB_ele";
    }
    else if (ele && !barrel) {
      name_InputHist = "hda_ele_EE"; name_InputHist_mllg = "hda_mllg_EE_ele";
      name_InputHist_ptllg = "hda_ptllg_EE_ele";
    }
    else if (!ele && barrel) {
      name_InputHist = "hda_mu_EB";name_InputHist_mllg = "hda_mllg_EB_mu";
      name_InputHist_ptllg = "hda_ptllg_EB_mu";
    }
    else if (!ele && !barrel) {
      name_InputHist = "hda_mu_EE"; name_InputHist_mllg = "hda_mllg_EE_mu";
      name_InputHist_ptllg = "hda_ptllg_EE_mu";
    }

    if (ele) name_InputHist_njet = "hda_njet_ele";
    else name_InputHist_njet = "hda_njet_mu";
  }
  else { //for mc
    if (ele && barrel) {
      name_InputHist = "hrec_ele_EB"; name_hgen ="hgen_ele_EB";
      name_InputHist_mllg = "hrec_mllg_EB_ele"; name_hgen_mllg = "hgen_mllg_EB_ele";
      name_InputHist_ptllg = "hrec_ptllg_EB_ele"; name_hgen_ptllg = "hgen_ptllg_EB_ele";}
    else if (ele && !barrel) {
      name_InputHist = "hrec_ele_EE"; name_hgen ="hgen_ele_EE";
      name_InputHist_mllg = "hrec_mllg_EE_ele"; name_hgen_mllg = "hgen_mllg_EE_ele";
      name_InputHist_ptllg = "hrec_ptllg_EB_ele"; name_hgen_ptllg = "hgen_ptllg_EB_ele";}
    else if (!ele && barrel) {
      name_InputHist = "hrec_mu_EB"; name_hgen ="hgen_mu_EB";
      name_InputHist_mllg = "hrec_mllg_EB_mu"; name_hgen_mllg = "hgen_mllg_EB_mu";
      name_InputHist_ptllg = "hrec_ptllg_EB_mu"; name_hgen_ptllg = "hgen_ptllg_EB_mu";}
    else  {
      name_InputHist = "hrec_mu_EE"; name_hgen ="hgen_mu_EE";
      name_InputHist_mllg = "hrec_mllg_EE_mu"; name_hgen_mllg = "hgen_mllg_EE_mu";
      name_InputHist_ptllg = "hrec_ptllg_EE_mu"; name_hgen_ptllg = "hgen_ptllg_EE_mu";}

    if (ele) {
      name_InputHist_njet = "hrec_njet_ele";
      name_hgen_njet = "hgen_njet_ele";
    }
    else {
      name_InputHist_njet = "hrec_njet_mu";
      name_hgen_njet = "hgen_njet_mu";}
  }

  TString name_hbkg_ttg_vv_ptg, name_hbkg_ttg_vv_mllg, name_hbkg_ttg_vv_ptllg, name_hbkg_ttg_vv_njet;
  TH1F *hbkg_ttg_vv_ptg, *hbkg_ttg_vv_mllg, *hbkg_ttg_vv_ptllg, *hbkg_ttg_vv_njet;

  //background
  if (data) {
    if (ele && barrel) {
      name_hbkg_ttg_vv_ptg = "hbkg_ttg_vv_ptg_EB_ele";
      name_hbkg_ttg_vv_mllg = "hbkg_ttg_vv_mllg_EB_ele";
      name_hbkg_ttg_vv_ptllg = "hbkg_ttg_vv_ptllg_EB_ele";
    }
    else if (ele && !barrel) {
      name_hbkg_ttg_vv_ptg = "hbkg_ttg_vv_ptg_EE_ele";
      name_hbkg_ttg_vv_mllg = "hbkg_ttg_vv_mllg_EE_ele";
      name_hbkg_ttg_vv_ptllg = "hbkg_ttg_vv_ptllg_EE_ele";
    }
    else if (!ele && barrel) {
      name_hbkg_ttg_vv_ptg = "hbkg_ttg_vv_ptg_EB_mu";
      name_hbkg_ttg_vv_mllg = "hbkg_ttg_vv_mllg_EB_mu";
      name_hbkg_ttg_vv_ptllg = "hbkg_ttg_vv_ptllg_EB_mu";
    }
    else if (!ele && !barrel) {
      name_hbkg_ttg_vv_ptg = "hbkg_ttg_vv_ptg_EE_mu";
      name_hbkg_ttg_vv_mllg = "hbkg_ttg_vv_mllg_EE_mu";
      name_hbkg_ttg_vv_ptllg = "hbkg_ttg_vv_ptllg_EE_mu";
    }
    if (ele) name_hbkg_ttg_vv_njet = "hbkg_ttg_vv_njet_ele";
    else     name_hbkg_ttg_vv_njet = "hbkg_ttg_vv_njet_mu";
  }

  histUnfoldInput = (TH1D*)finput->Get(name_InputHist);
  histUnfoldInput_mllg = (TH1D*) finput->Get(name_InputHist_mllg);
  histUnfoldInput_ptllg = (TH1D*) finput->Get(name_InputHist_ptllg);
  histUnfoldInput_njet = (TH1D*) finput->Get(name_InputHist_njet);

  hbkg_ttg_vv_ptg = (TH1F*) finput->Get(name_hbkg_ttg_vv_ptg);
  hbkg_ttg_vv_mllg = (TH1F*) finput->Get(name_hbkg_ttg_vv_mllg);
  hbkg_ttg_vv_ptllg = (TH1F*) finput->Get(name_hbkg_ttg_vv_ptllg);
  hbkg_ttg_vv_njet = (TH1F*) finput->Get(name_hbkg_ttg_vv_njet);

  if (!data) {
    hgen = (TH1D*) finput->Get(name_hgen);
    hgen_mllg = (TH1D*) finput->Get(name_hgen_mllg);
    hgen_njet = (TH1D*) finput->Get(name_hgen_njet);
  }


  cout << "done reading input" << endl;

  histUnfoldInput->Sumw2();
  histUnfoldInput_mllg->Sumw2();
  histUnfoldInput_ptllg->Sumw2();
  histUnfoldInput_njet->Sumw2();

  if (data) {
    hbkg_ttg_vv_ptg->Sumw2();
    hbkg_ttg_vv_mllg->Sumw2();
    hbkg_ttg_vv_ptllg->Sumw2();
    hbkg_ttg_vv_njet->Sumw2();
  }

  cout << "------- unfolding photon pt -------------" << endl;
  const Double_t* xval_Det = histUnfoldInput->GetXaxis()->GetXbins()->GetArray();
  const Double_t* xval_Gen = histUnfoldMatrix->GetXaxis()->GetXbins()->GetArray();

  int nbins    = histUnfoldInput ->GetXaxis()->GetNbins();
  int ngenbins = histUnfoldMatrix->GetXaxis()->GetNbins();

  cout << " reco bins " << nbins  << endl;
  cout << "gen bins " << ngenbins << endl;


  for (int i = 0; i < nbins; ++i)
    {
      cout << xval_Det[i] << "\t";
    }
  cout<< xval_Det[nbins] << endl;

  double genpt[ngenbins+1];

  for (int i = 0; i < ngenbins; ++i)
    {
      cout << xval_Gen[i] << "\t";
    }

  cout<< xval_Gen[ngenbins] << endl;
  TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature; //kRegModeDerivative;
  TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
  TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeNone; //TUnfoldDensity::kDensityModeBinWidth;

  TUnfoldDensity unfold(histUnfoldMatrix,TUnfold::kHistMapOutputHoriz, regMode,constraintMode, densityFlags);

  if(unfold.SetInput(histUnfoldInput)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }

  if (data) {
    unfold.SubtractBackground(hbkg_ttg_vv_ptg, "ttg_vv_bkg", 1, 0);
  }
  cout << "running the unfolding" << endl;

  //Int_t nScan=30; 
  Double_t tauMin=0;
  Double_t tauMax=0;

  Int_t nScan=50;
  //Double_t tauMin=1.E-10;
  //Double_t tauMax=0.01;
  TSpline *logTauX,*logTauY, *logTauCurvature;
  TGraph *lCurve;

  Int_t iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);

  cout << "tau = " << unfold.GetTau() << endl;
  cout<<"chi**2 = "<<unfold.GetChi2A()<<" t Chi2L: "<<unfold.GetChi2L()
      <<" / "<<unfold.GetNdf()<<"\n";

  Double_t t[1],x[1],y[1];
  logTauX->GetKnot(iBest,t[0],x[0]);
  logTauY->GetKnot(iBest,t[0],y[0]);
  TGraph *bestLcurve=new TGraph(1,x,y);
  TGraph *bestLogTauLogChi2=new TGraph(1,t,x);


  TH1 *histUnfoldOutput= unfold.GetOutput("histUnfoldOutput");
  TH1 *histFolded=unfold.GetFoldedOutput("histFolded");

  TH2 *histEmatStat=unfold.GetEmatrixInput("histEmatStat");
  TH2 *histEmatTotal=unfold.GetEmatrixTotal("histEmatTotal");
  TH1 *hbias = unfold.GetBias("hbias");
  TH1 *hinput = unfold.GetInput("hinput");
  TH2 *hProb = unfold.GetProbabilityMatrix("hProb");
  TH2 *hmatrixUnCorr = unfold.GetEmatrixSysUncorr("hmatrixUnCorr");

  TH1 *histUnfoldStat=new TH1D("PT(unfold,staterr)",";Pt(gen)",ngenbins,xval_Gen);
  TH1 *histUnfoldTotal=new TH1D("PT(unfold,totalerr)",";Pt(gen)",ngenbins,xval_Gen);

  cout << "looooooooooooooooooooooooooooooop " << endl;
  for(Int_t i=0;i<ngenbins+2;i++) {
    Double_t c=histUnfoldOutput->GetBinContent(i);

    // histogram with unfolded data and stat errors
    histUnfoldStat->SetBinContent(i,c);
    histUnfoldStat->SetBinError(i,TMath::Sqrt(histEmatStat->GetBinContent(i,i)));

    // histogram with unfolded data and total errors
    histUnfoldTotal->SetBinContent(i,c);
    histUnfoldTotal->SetBinError(i,TMath::Sqrt(histEmatTotal->GetBinContent(i,i)));
  }

  TH2 *gHistInvEMatrix;
  // get global correlation coefficients
  TH1 *histRhoi=unfold.GetRhoItotal("rho_I",
                                    0, // use default title
                                    0, // all distributions
                                    "", // discard underflow and overflow bins on all axes
                                    kTRUE, // use original binning
                                    &gHistInvEMatrix // store inverse of error matrix
                                    );


  TH2 *histCorr = unfold.GetRhoIJtotal("histCorr");

  //================================================================================//
  // unfolding Mllg
  const Double_t* xval_Det_mllg = histUnfoldInput_mllg->GetXaxis()->GetXbins()->GetArray();
  const Double_t* xval_Gen_mllg = histUnfoldMatrix_mllg->GetXaxis()->GetXbins()->GetArray();

  int nbins_mllg    = histUnfoldInput_mllg ->GetXaxis()->GetNbins();
  int ngenbins_mllg = histUnfoldMatrix_mllg->GetXaxis()->GetNbins();

  cout << "---------binning of mllg-----------" << endl;
  for (int i = 0; i < nbins_mllg; ++i)
    {
      cout << xval_Det_mllg[i] << "\t";
    }
  cout<< xval_Det_mllg[nbins_mllg] << endl;

  double genmllg[ngenbins_mllg+1];

  for (int i = 0; i < ngenbins_mllg; ++i)
    {
      cout << xval_Gen_mllg[i] << "\t";
    }
  cout<< xval_Gen_mllg[ngenbins_mllg] << endl;


  TUnfoldDensity unfold_mllg(histUnfoldMatrix_mllg,TUnfold::kHistMapOutputHoriz, regMode,constraintMode, densityFlags);

  if(unfold_mllg.SetInput(histUnfoldInput_mllg)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }

  if (data) {
    unfold_mllg.SubtractBackground(hbkg_ttg_vv_mllg, "ttg_vv_bkg", 1, 0);
  }

  TSpline *logTauX_mllg,*logTauY_mllg, *logTauCurvature_mllg;
  TGraph *lCurve_mllg;
  TSpline *rhoLogTau=0;
  const char *SCAN_DISTRIBUTION="signal";
  const char *SCAN_AXISSTEERING=0;
  Int_t iBest_mllg=unfold_mllg.ScanLcurve(nScan,tauMin,tauMax,&lCurve_mllg,&logTauX_mllg,&logTauY_mllg,&logTauCurvature_mllg);

  cout << "tau = " << unfold_mllg.GetTau() << endl;
  cout<<"chi**2 = "<<unfold_mllg.GetChi2A()<<" t Chi2L: "<<unfold_mllg.GetChi2L()
      <<" / "<<unfold_mllg.GetNdf()<<"\n";

  TH1 *histUnfoldOutput_mllg= unfold_mllg.GetOutput("histUnfoldOutput_mllg");
  TH1 *histFolded_mllg=unfold_mllg.GetFoldedOutput("histFolded_mllg");
  TH2 *histEmatStat_mllg=unfold_mllg.GetEmatrixInput("histEmatStat_mllg");
  TH2 *histEmatTotal_mllg=unfold_mllg.GetEmatrixTotal("histEmatTotal_mllg");
  TH1 *hbias_mllg = unfold_mllg.GetBias("hbias_mllg");
  TH1 *hinput_mllg = unfold_mllg.GetInput("hinput_mllg");
  TH2 *hProb_mllg = unfold_mllg.GetProbabilityMatrix("hProb_mllg");
  TH2 *hmatrixUnCorr_mllg = unfold_mllg.GetEmatrixSysUncorr("hmatrixUnCorr_mllg");

  TH1 *histUnfoldStat_mllg=new TH1D("Mllg(unfold_mllg,staterr)",";Mllg(gen)",ngenbins_mllg,xval_Gen_mllg);
  TH1 *histUnfoldTotal_mllg=new TH1D("Mllg(unfold_mllg,totalerr)",";Mllg(gen)",ngenbins_mllg,xval_Gen_mllg);

  for(Int_t i=0;i<ngenbins_mllg+2;i++) {
    Double_t c=histUnfoldOutput_mllg->GetBinContent(i);

    histUnfoldStat_mllg->SetBinContent(i,c);
    histUnfoldStat_mllg->SetBinError(i,TMath::Sqrt(histEmatStat_mllg->GetBinContent(i,i)));

    histUnfoldTotal_mllg->SetBinContent(i,c);
    histUnfoldTotal_mllg->SetBinError(i,TMath::Sqrt(histEmatTotal_mllg->GetBinContent(i,i)));
  }

  TH2 *gHistInvEMatrix_mllg;
  TH1 *histRhoi_mllg = unfold_mllg.GetRhoItotal("rho_I_mllg", 0, 0, "", kTRUE, &gHistInvEMatrix_mllg);

  TH2 *histCorr_mllg = unfold_mllg.GetRhoIJtotal("histCorr_mllg");


  //================================================================================//
  // unfolding ptllg
  const Double_t* xval_Det_ptllg = histUnfoldInput_ptllg->GetXaxis()->GetXbins()->GetArray();
  const Double_t* xval_Gen_ptllg = histUnfoldMatrix_ptllg->GetXaxis()->GetXbins()->GetArray();

  int nbins_ptllg    = histUnfoldInput_ptllg ->GetXaxis()->GetNbins();
  int ngenbins_ptllg = histUnfoldMatrix_ptllg->GetXaxis()->GetNbins();

  cout << "---------binning of ptllg-----------" << endl;
  for (int i = 0; i < nbins_ptllg; ++i)
    {
      cout << xval_Det_ptllg[i] << "\t";
    }
  cout<< xval_Det_ptllg[nbins_ptllg] << endl;

  double genptllg[ngenbins_ptllg+1];

  for (int i = 0; i < ngenbins_ptllg; ++i)
    {
      cout << xval_Gen_ptllg[i] << "\t";
    }
  cout<< xval_Gen_ptllg[ngenbins_ptllg] << endl;


  TUnfoldDensity unfold_ptllg(histUnfoldMatrix_ptllg,TUnfold::kHistMapOutputHoriz, regMode,constraintMode, densityFlags);

  if(unfold_ptllg.SetInput(histUnfoldInput_ptllg)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }

  //if (data) {
  // unfold_ptllg.SubtractBackground(hbkg_ttg_vv_ptllg, "ttg_vv_bkg", 1, 0);
  //}

  TSpline *logTauX_ptllg,*logTauY_ptllg, *logTauCurvature_ptllg;
  TGraph *lCurve_ptllg;
  Int_t iBest_ptllg=unfold_ptllg.ScanLcurve(nScan,tauMin,tauMax,&lCurve_ptllg,&logTauX_ptllg,&logTauY_ptllg,&logTauCurvature_ptllg);

  cout << "tau = " << unfold_ptllg.GetTau() << endl;
  cout<<"chi**2 = "<<unfold_ptllg.GetChi2A()<<" t Chi2L: "<<unfold_ptllg.GetChi2L()
      <<" / "<<unfold_ptllg.GetNdf()<<"\n";

  TH1 *histUnfoldOutput_ptllg= unfold_ptllg.GetOutput("histUnfoldOutput_ptllg");
  TH1 *histFolded_ptllg=unfold_ptllg.GetFoldedOutput("histFolded_ptllg");
  TH2 *histEmatStat_ptllg=unfold_ptllg.GetEmatrixInput("histEmatStat_ptllg");
  TH2 *histEmatTotal_ptllg=unfold_ptllg.GetEmatrixTotal("histEmatTotal_ptllg");
  TH1 *hbias_ptllg = unfold_ptllg.GetBias("hbias_ptllg");
  TH1 *hinput_ptllg = unfold_ptllg.GetInput("hinput_ptllg");
  TH2 *hProb_ptllg = unfold_ptllg.GetProbabilityMatrix("hProb_ptllg");
  TH2 *hmatrixUnCorr_ptllg = unfold_ptllg.GetEmatrixSysUncorr("hmatrixUnCorr_ptllg");

  TH1 *histUnfoldStat_ptllg=new TH1D("Ptllg(unfold_ptllg,staterr)",";Ptllg(gen)",ngenbins_ptllg,xval_Gen_ptllg);
  TH1 *histUnfoldTotal_ptllg=new TH1D("Ptllg(unfold_ptllg,totalerr)",";Ptllg(gen)",ngenbins_ptllg,xval_Gen_ptllg);

  for(Int_t i=0;i<ngenbins_ptllg+2;i++) {
    Double_t c=histUnfoldOutput_ptllg->GetBinContent(i);

    histUnfoldStat_ptllg->SetBinContent(i,c);
    histUnfoldStat_ptllg->SetBinError(i,TMath::Sqrt(histEmatStat_ptllg->GetBinContent(i,i)));

    histUnfoldTotal_ptllg->SetBinContent(i,c);
    histUnfoldTotal_ptllg->SetBinError(i,TMath::Sqrt(histEmatTotal_ptllg->GetBinContent(i,i)));
  }

  TH2 *gHistInvEMatrix_ptllg;
  TH1 *histRhoi_ptllg = unfold_ptllg.GetRhoItotal("rho_I_ptllg", 0, 0, "", kTRUE, &gHistInvEMatrix_ptllg);

  TH2 *histCorr_ptllg = unfold_ptllg.GetRhoIJtotal("histCorr_ptllg");


  //=====================================//
  // jet multiplicity
  //const Double_t* xval_Det_njet = histUnfoldInput_njet->GetXaxis()->GetXbins()->GetArray();

  int nbins_njet    = histUnfoldInput_njet ->GetXaxis()->GetNbins();
  int ngenbins_njet = histUnfoldMatrix_njet->GetXaxis()->GetNbins();

  float xval_Gen_njet[ngenbins_njet+1];
  float xval_Det_njet[nbins_njet+1];

  for (int i = 0; i <= nbins_njet; ++i)
    {
      if (i == nbins_njet) xval_Det_njet[i] = histUnfoldMatrix_njet->GetYaxis()->GetBinUpEdge(i);
      else  xval_Det_njet[i] = histUnfoldMatrix_njet->GetYaxis()->GetBinLowEdge(i+1);
      cout << xval_Det_njet[i] << "\t";

    }
  cout << endl;


  cout << "nbins of njet: " << ngenbins_njet << endl;
  for (int i = 0; i <= ngenbins_njet; ++i)
    {
      if (i == ngenbins_njet) xval_Gen_njet[i] = histUnfoldMatrix_njet->GetXaxis()->GetBinUpEdge(i);
      else  xval_Gen_njet[i] = histUnfoldMatrix_njet->GetXaxis()->GetBinLowEdge(i+1);
      cout << xval_Gen_njet[i] << "\t";
    }
  cout<< endl;

  TUnfoldDensity unfold_njet(histUnfoldMatrix_njet,TUnfold::kHistMapOutputHoriz, regMode,constraintMode, densityFlags);

  if(unfold_njet.SetInput(histUnfoldInput_njet)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }

  if (data) {
    unfold_njet.SubtractBackground(hbkg_ttg_vv_njet, "ttg_vv_bkg", 1, 0);
  }

  TSpline *logTauX_njet,*logTauY_njet, *logTauCurvature_njet;
  TGraph *lCurve_njet;
  //Int_t iBest_njet=unfold_njet.ScanLcurve(nScan,tauMin,tauMax,&lCurve_njet,&logTauX_njet,&logTauY_njet,&logTauCurvature_njet);
  Int_t iBest_njet=unfold_njet.ScanLcurve(nScan,tauMin,tauMax,&lCurve_njet,&logTauX_njet,&logTauY_njet);

  TH1 *histUnfoldOutput_njet= unfold_njet.GetOutput("histUnfoldOutput_njet");
  TH1 *histFolded_njet=unfold_njet.GetFoldedOutput("histFolded_njet");
  TH2 *histEmatStat_njet=unfold_njet.GetEmatrixInput("histEmatStat_njet");
  TH2 *histEmatTotal_njet=unfold_njet.GetEmatrixTotal("histEmatTotal_njet");
  TH1 *hbias_njet = unfold_njet.GetBias("hbias_njet");
  TH1 *hinput_njet = unfold_njet.GetInput("hinput_njet");
  TH2 *hProb_njet = unfold_njet.GetProbabilityMatrix("hProb_njet");
  TH2 *hmatrixUnCorr_njet = unfold_njet.GetEmatrixSysUncorr("hmatrixUnCorr_njet");

  TH1 *histUnfoldStat_njet=new TH1D("Njet(unfold_njet,staterr)",";Njet(gen)",ngenbins_njet,xval_Gen_njet);
  TH1 *histUnfoldTotal_njet=new TH1D("Njet(unfold_njet,totalerr)",";Njet(gen)",ngenbins_njet,xval_Gen_njet);

  for(Int_t i=0;i<ngenbins_njet+2;i++) {
    Double_t c=histUnfoldOutput_njet->GetBinContent(i);

    histUnfoldStat_njet->SetBinContent(i,c);
    histUnfoldStat_njet->SetBinError(i,TMath::Sqrt(histEmatStat_njet->GetBinContent(i,i)));

    histUnfoldTotal_njet->SetBinContent(i,c);
    histUnfoldTotal_njet->SetBinError(i,TMath::Sqrt(histEmatTotal_njet->GetBinContent(i,i)));
  }

  TH2 *gHistInvEMatrix_njet;
  TH1 *histRhoi_njet = unfold_njet.GetRhoItotal("rho_I_njet", 0, 0, "", kTRUE, &gHistInvEMatrix_njet);

  TH2 *histCorr_njet = unfold_njet.GetRhoIJtotal("histCorr_njet");

  TString filename="outputUnfold/SB_EE7to13_EE6to14/Zg_Unfolding_Subtract_Bkg_";
  if (data) filename += "data_";
  else filename +="MC_";
  if (ele) filename +="ele_";
  else filename +="mu_";
  if (barrel) filename +="barrel_";
  else filename +="endcap_";
  //filename += "biasDn_";
  //filename += "ttg_vv_bkg_UncDn_";
  filename += suf;

  TFile* fout = new TFile(filename,"RECREATE");

  histUnfoldInput->Write("histUnfoldInput");
  histUnfoldStat->Write("histUnfoldStat");
  histUnfoldOutput->Write("histUnfoldOutput");
  histUnfoldTotal->Write("histUnfoldTotal");
  histUnfoldMatrix->Write("histUnfoldMatrix");
  bestLogTauLogChi2->Write("bestLogTauLogChi2");
  bestLcurve->Write("bestLcurve");
  hbias->Write("hbias");
  hProb->Write("hProb");
  hmatrixUnCorr->Write("hmatrixUnCorr");
  histCorr->Write("histCorr");
  histFolded->Write("histFolded");
  //lCurve->Write("lCurve");
  //logTauX->Write("logTauX");
  //logTauY->Write("logTauY");
  //logTauCurvature->Write("logTauCurvature"); 

  histUnfoldInput_mllg->Write("histUnfoldInput_mllg");
  histUnfoldOutput_mllg->Write("histUnfoldOutput_mllg");
  histUnfoldStat_mllg->Write("histUnfoldStat_mllg");
  histUnfoldTotal_mllg->Write("histUnfoldTotal_mllg");
  histUnfoldMatrix_mllg->Write("histUnfoldMatrix_mllg");
  hProb_mllg->Write("hProb_mllg");
  histCorr_mllg->Write("histCorr_mllg");
  histFolded_mllg->Write("histFolded_mllg");
  lCurve_mllg->Write("lCurve_mllg");
  logTauX_mllg->Write("logTauX_mllg");
  logTauY_mllg->Write("logTauY_mllg");
  logTauCurvature_mllg->Write("logTauCurvature_mllg");

  histUnfoldInput_ptllg->Write("histUnfoldInput_ptllg");
  histUnfoldOutput_ptllg->Write("histUnfoldOutput_ptllg");
  histUnfoldStat_ptllg->Write("histUnfoldStat_ptllg");
  histUnfoldTotal_ptllg->Write("histUnfoldTotal_ptllg");
  histUnfoldMatrix_ptllg->Write("histUnfoldMatrix_ptllg");
  hProb_ptllg->Write("hProb_ptllg");
  histCorr_ptllg->Write("histCorr_ptllg");
  histFolded_ptllg->Write("histFolded_ptllg");
  lCurve_ptllg->Write("lCurve_ptllg");
  logTauX_ptllg->Write("logTauX_ptllg");
  logTauY_ptllg->Write("logTauY_ptllg");
  logTauCurvature_ptllg->Write("logTauCurvature_ptllg");


  histUnfoldInput_njet->Write("histUnfoldInput_njet");
  histUnfoldOutput_njet->Write("histUnfoldOutput_njet");
  histUnfoldStat_njet->Write("histUnfoldStat_njet");
  histUnfoldTotal_njet->Write("histUnfoldTotal_njet");
  histUnfoldMatrix_njet->Write("histUnfoldMatrix_njet");
  hProb_njet->Write("hProb_njet");
  histCorr_njet->Write("histCorr_njet");
  histFolded_njet->Write("histFolded_njet");
  //lCurve_njet->Write("lCurve_njet");
  //logTauX_njet->Write("logTauX_njet");
  //logTauY_njet->Write("logTauY_njet"); 
  //logTauCurvature_njet->Write("logTauCurvature_njet");

  if (!data) {
    hgen->Write("hgen");
    hgen_mllg->Write("hgen_mllg");
    hgen_ptllg->Write("hgen_ptllg");
    hgen_njet->Write("hgen_njet");
  }

}

