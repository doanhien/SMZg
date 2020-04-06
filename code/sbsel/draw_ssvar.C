#include "TFile.h"
#include "TH1.h"

void draw_ssvar() {


  gStyle->SetOptStat(0);

  TFile *fdy = new TFile("../ana/minitrees/ZJets_aMCatNLO_Summer16_TMVA420_UpTp6000_5VarCorr_etawei.root", "read");
  TFile *fqcd = new TFile("input/bkgTemplate_summer16_420_BDT_Up6000_5VarCorr_job_summer16_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root", "read");
  TFile *fda = new TFile("../ana/minitrees/DoubleEG_Run2016_FebReminiAOD_Summer16_TMVA420_UpTp6000_5VarCorr_etawei.root", "read");

  TTree *tdy = (TTree*) fdy->Get("outtree");
  TTree *tda = (TTree*) fda->Get("outtree");
  TTree *tqcd = (TTree*) fqcd->Get("outtree");

  TH1F *hetawidth_dy = new TH1F("hetawidth_dy", "", 25, 0., 0.1);
  TH1F *hphiwidth_dy = new TH1F("hphiwidth_dy", "", 40, 0., 0.2);
  TH1F *hscRawE_dy   = new TH1F("hscRawE_dy", "", 50, 0, 200);
  TH1F *hs4_dy = new TH1F("hs4_dy", "", 20, 0., 1.2);
  TH1F *hr9_dy = new TH1F("hr9_dy", "", 20, 0., 1.2);
  TH1F *hsieip_dy = new TH1F("hsieip_dy", "", 50, -0.001, 0.001);
  TH1F *hsieie_dy = new TH1F("hsieie_dy", "", 50, 0, 0.05);
  TH1F *hesrr_dy = new TH1F("hesrr_dy", "", 50, 0, 15);
  TH1F *hEsOEn_dy = new TH1F("hEsOEn_dy", "", 40, 0, 0.8);
  TH1F *heta_dy = new TH1F("heta_dy", "", 25, -2.5, 2.5);

  TH1F *hetawidth_dy_sb = new TH1F("hetawidth_dy_sb", "", 25, 0., 0.1);
  TH1F *hphiwidth_dy_sb = new TH1F("hphiwidth_dy_sb", "", 40, 0., 0.2);
  TH1F *hscRawE_dy_sb   = new TH1F("hscRawE_dy_sb", "", 50, 0, 200);
  TH1F *hs4_dy_sb = new TH1F("hs4_dy_sb", "", 20, 0., 1.2);
  TH1F *hr9_dy_sb = new TH1F("hr9_dy_sb", "", 20, 0., 1.2);
  TH1F *hsieip_dy_sb = new TH1F("hsieip_dy_sb", "", 50, -0.001, 0.001);
  TH1F *hsieie_dy_sb = new TH1F("hsieie_dy_sb", "", 50, 0, 0.05);
  TH1F *hesrr_dy_sb = new TH1F("hesrr_dy_sb", "", 50, 0, 15);
  TH1F *hEsOEn_dy_sb = new TH1F("hEsOEn_dy_sb", "", 40, 0, 0.8);
  TH1F *heta_dy_sb = new TH1F("heta_dy_sb", "", 25, -2.5, 2.5);

  TH1F *hetawidth_da_sb = new TH1F("hetawidth_da_sb", "", 25, 0., 0.1);
  TH1F *hphiwidth_da_sb = new TH1F("hphiwidth_da_sb", "", 40, 0., 0.2);
  TH1F *hscRawE_da_sb   = new TH1F("hscRawE_da_sb", "", 50, 0, 200);
  TH1F *hs4_da_sb = new TH1F("hs4_da_sb", "", 20, 0., 1.2);
  TH1F *hr9_da_sb = new TH1F("hr9_da_sb", "", 20, 0., 1.2);
  TH1F *hsieip_da_sb = new TH1F("hsieip_da_sb", "", 50, -0.001, 0.001);
  TH1F *hsieie_da_sb = new TH1F("hsieie_da_sb", "", 50, 0, 0.05);
  TH1F *hesrr_da_sb = new TH1F("hesrr_da_sb", "", 50, 0, 15);
  TH1F *hEsOEn_da_sb = new TH1F("hEsOEn_da_sb", "", 40, 0, 0.8);
  TH1F *heta_da_sb = new TH1F("heta_da_sb", "", 25, -2.5, 2.5);


  TH1F *hetawidth_qcd_sb = new TH1F("hetawidth_qcd_sb", "", 25, 0., 0.1);
  TH1F *hphiwidth_qcd_sb = new TH1F("hphiwidth_qcd_sb", "", 40, 0., 0.2);
  TH1F *hscRawE_qcd_sb   = new TH1F("hscRawE_qcd_sb", "", 50, 0, 200);
  TH1F *hs4_qcd_sb = new TH1F("hs4_qcd_sb", "", 20, 0., 1.2);
  TH1F *hr9_qcd_sb = new TH1F("hr9_qcd_sb", "", 20, 0., 1.2);
  TH1F *hsieip_qcd_sb = new TH1F("hsieip_qcd_sb", "", 50, -0.001, 0.001);
  TH1F *hsieie_qcd_sb = new TH1F("hsieie_qcd_sb", "", 50, 0, 0.05);
  TH1F *hesrr_qcd_sb = new TH1F("hesrr_qcd_sb", "", 50, 0, 15);
  TH1F *hEsOEn_qcd_sb = new TH1F("hEsOEn_qcd_sb", "", 40, 0, 0.8);
  TH1F *heta_qcd_sb = new TH1F("heta_qcd_sb", "", 25, -2.5, 2.5);
  

  TCut cut = "isEE==1 && gamma_ChIso<1.5 && gamma_pt>15 && gamma_pt<20 && leptType==11 && trig_Ele23_Ele12==1";
  TCut cut_sb = "isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso < 13. && ((leptType==11 && trig_Ele23_Ele12==1) || (leptType==13 && trig_Mu17_Mu8==1))";
  TCut cut_sb_qcd = "isEE==1 && gamma_pt>15 && gamma_pt<20 && gamma_ChIso>7. && gamma_ChIso < 13.";
  TCut wei = "puweigj_65nb*genWeight";
  TCut wei_qcd = "puwei*genwei";

  //signal region
  tdy->Draw("gamma_SCEtaWidth_rw >> hetawidth_dy", cut*wei);
  tdy->Draw("gamma_scphiwidth >> hphiwidth_dy", cut*wei);
  tdy->Draw("gamma_SCRawE >> hscRawE_dy", cut*wei);
  tdy->Draw("gamma_s4Full5x5_rw >> hs4_dy", cut*wei);
  tdy->Draw("gamma_R9_rw >> hr9_dy", cut*wei);
  tdy->Draw("gamma_sigmaIetaIphi_rw >> hsieip_dy", cut*wei);
  tdy->Draw("gamma_sigmaIetaIeta_rw >> hsieie_dy", cut*wei);
  tdy->Draw("gamma_esRR >> hesrr_dy", cut*wei);
  tdy->Draw("(gamma_ESEnP1+gamma_ESEnP2)/gamma_SCRawE >> hEsOEn_dy", cut*wei);
  tdy->Draw("gamma_sceta >> heta_dy", cut*wei);

  //side-band region
  tdy->Draw("gamma_SCEtaWidth_rw >> hetawidth_dy_sb", cut_sb*wei*"etawei");
  tdy->Draw("gamma_scphiwidth >> hphiwidth_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_SCRawE >> hscRawE_dy_sb", cut_sb*wei*"etawei");
  tdy->Draw("gamma_s4Full5x5_rw >> hs4_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_R9_rw >> hr9_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_sigmaIetaIphi_rw >> hsieip_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_sigmaIetaIeta_rw >> hsieie_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_esRR >> hesrr_dy_sb", cut_sb*wei);
  tdy->Draw("(gamma_ESEnP1+gamma_ESEnP2)/gamma_SCRawE >> hEsOEn_dy_sb", cut_sb*wei);
  tdy->Draw("gamma_sceta >> heta_dy_sb", cut_sb*wei*"etawei");

  //data side-band
  tda->Draw("gamma_SCEtaWidth >> hetawidth_da_sb", cut_sb*"etawei");
  tda->Draw("gamma_scphiwidth >> hphiwidth_da_sb", cut_sb);
  tda->Draw("gamma_SCRawE >> hscRawE_da_sb", cut_sb*"etawei");
  tda->Draw("gamma_s4Full5x5 >> hs4_da_sb", cut_sb);
  tda->Draw("gamma_R9 >> hr9_da_sb", cut_sb);
  tda->Draw("gamma_sigmaIetaIphi >> hsieip_da_sb", cut_sb);
  tda->Draw("gamma_sigmaIetaIeta >> hsieie_da_sb", cut_sb);
  tda->Draw("gamma_esRR >> hesrr_da_sb", cut_sb);
  tda->Draw("(gamma_ESEnP1+gamma_ESEnP2)/gamma_SCRawE >> hEsOEn_da_sb", cut_sb);
  tda->Draw("gamma_sceta >> heta_da_sb", cut_sb*"etawei");

  //qcd side-band
  tqcd->Draw("gamma_SCEtaWidth_rw >> hetawidth_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_scphiwidth >> hphiwidth_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_SCRawE >> hscRawE_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_s4Full5x5_rw >> hs4_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_R9_rw >> hr9_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_sigmaIetaIphi_rw >> hsieip_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_sigmaIetaIeta_rw >> hsieie_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_esRR >> hesrr_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("(gamma_ESEnP1+gamma_ESEnP2)/gamma_SCRawE >> hEsOEn_qcd_sb", cut_sb_qcd*wei_qcd);
  tqcd->Draw("gamma_eta >> heta_qcd_sb", cut_sb_qcd*wei_qcd);

  //clone histogram for error display
  hetawidth_dy_err = (TH1F*) hetawidth_dy->Clone();
  hphiwidth_dy_err = (TH1F*) hphiwidth_dy->Clone();
  hscRawE_dy_err = (TH1F*) hscRawE_dy->Clone();
  hs4_dy_err = (TH1F*) hs4_dy->Clone();
  hr9_dy_err = (TH1F*) hr9_dy->Clone();
  hsieip_dy_err = (TH1F*) hsieip_dy->Clone();
  hsieie_dy_err = (TH1F*) hsieie_dy->Clone();
  hesrr_dy_err = (TH1F*) hesrr_dy->Clone();
  hEsOEn_dy_err = (TH1F*) hEsOEn_dy->Clone();
  heta_dy_err = (TH1F*) heta_dy->Clone();

  hetawidth_dy_sb_err = (TH1F*) hetawidth_dy_sb->Clone();
  hphiwidth_dy_sb_err = (TH1F*) hphiwidth_dy_sb->Clone();
  hscRawE_dy_sb_err = (TH1F*) hscRawE_dy_sb->Clone();
  hs4_dy_sb_err = (TH1F*) hs4_dy_sb->Clone();
  hr9_dy_sb_err = (TH1F*) hr9_dy_sb->Clone();
  hsieip_dy_sb_err = (TH1F*) hsieip_dy_sb->Clone();
  hsieie_dy_sb_err = (TH1F*) hsieie_dy_sb->Clone();
  hesrr_dy_sb_err = (TH1F*) hesrr_dy_sb->Clone();
  hEsOEn_dy_sb_err = (TH1F*) hEsOEn_dy_sb->Clone();
  heta_dy_sb_err = (TH1F*) heta_dy_sb->Clone();

  hetawidth_qcd_sb_err = (TH1F*) hetawidth_qcd_sb->Clone();
  hphiwidth_qcd_sb_err = (TH1F*) hphiwidth_qcd_sb->Clone();
  hscRawE_qcd_sb_err = (TH1F*) hscRawE_qcd_sb->Clone();
  hs4_qcd_sb_err = (TH1F*) hs4_qcd_sb->Clone();
  hr9_qcd_sb_err = (TH1F*) hr9_qcd_sb->Clone();
  hsieip_qcd_sb_err = (TH1F*) hsieip_qcd_sb->Clone();
  hsieie_qcd_sb_err = (TH1F*) hsieie_qcd_sb->Clone();
  hesrr_qcd_sb_err = (TH1F*) hesrr_qcd_sb->Clone();
  hEsOEn_qcd_sb_err = (TH1F*) hEsOEn_qcd_sb->Clone();
  heta_qcd_sb_err = (TH1F*) heta_qcd_sb->Clone();


  //normalize to 1
  hetawidth_dy_err->Scale(1./hetawidth_dy_err->Integral());
  hphiwidth_dy_err->Scale(1./hphiwidth_dy_err->Integral());
  hscRawE_dy_err->Scale(1./hscRawE_dy_err->Integral());
  hs4_dy_err->Scale(1./hs4_dy_err->Integral());
  hr9_dy_err->Scale(1./hr9_dy_err->Integral());
  hsieip_dy_err->Scale(1./hsieip_dy_err->Integral());
  hsieie_dy_err->Scale(1./hsieie_dy_err->Integral());
  hesrr_dy_err->Scale(1./hesrr_dy_err->Integral());
  hEsOEn_dy_err->Scale(1./hEsOEn_dy_err->Integral());
  heta_dy_err->Scale(1./heta_dy_err->Integral());

  hetawidth_dy_sb_err->Scale(1./hetawidth_dy_sb_err->Integral());
  hphiwidth_dy_sb_err->Scale(1./hphiwidth_dy_sb_err->Integral());
  hscRawE_dy_sb_err->Scale(1./hscRawE_dy_sb_err->Integral());
  hs4_dy_sb_err->Scale(1./hs4_dy_sb_err->Integral());
  hr9_dy_sb_err->Scale(1./hr9_dy_sb_err->Integral());
  hsieip_dy_sb_err->Scale(1./hsieip_dy_sb_err->Integral());
  hsieie_dy_sb_err->Scale(1./hsieie_dy_sb_err->Integral());
  hesrr_dy_sb_err->Scale(1./hesrr_dy_sb_err->Integral());
  hEsOEn_dy_sb_err->Scale(1./hEsOEn_dy_sb_err->Integral());
  heta_dy_sb_err->Scale(1./heta_dy_sb_err->Integral());

  hetawidth_qcd_sb_err->Scale(1./hetawidth_qcd_sb_err->Integral());
  hphiwidth_qcd_sb_err->Scale(1./hphiwidth_qcd_sb_err->Integral());
  hscRawE_qcd_sb_err->Scale(1./hscRawE_qcd_sb_err->Integral());
  hs4_qcd_sb_err->Scale(1./hs4_qcd_sb_err->Integral());
  hr9_qcd_sb_err->Scale(1./hr9_qcd_sb_err->Integral());
  hsieip_qcd_sb_err->Scale(1./hsieip_qcd_sb_err->Integral());
  hsieie_qcd_sb_err->Scale(1./hsieie_qcd_sb_err->Integral());
  hesrr_qcd_sb_err->Scale(1./hesrr_qcd_sb_err->Integral());
  hEsOEn_qcd_sb_err->Scale(1./hEsOEn_qcd_sb_err->Integral());
  heta_qcd_sb_err->Scale(1./heta_qcd_sb_err->Integral());


  //normalize to 1
  hetawidth_dy->Scale(1./hetawidth_dy->Integral());
  hphiwidth_dy->Scale(1./hphiwidth_dy->Integral());
  hscRawE_dy->Scale(1./hscRawE_dy->Integral());
  hs4_dy->Scale(1./hs4_dy->Integral());
  hr9_dy->Scale(1./hr9_dy->Integral());
  hsieip_dy->Scale(1./hsieip_dy->Integral());
  hsieie_dy->Scale(1./hsieie_dy->Integral());
  hesrr_dy->Scale(1./hesrr_dy->Integral());
  hEsOEn_dy->Scale(1./hEsOEn_dy->Integral());
  heta_dy->Scale(1./heta_dy->Integral());

  
  hetawidth_dy_sb->Scale(1./hetawidth_dy_sb->Integral());
  hphiwidth_dy_sb->Scale(1./hphiwidth_dy_sb->Integral());
  hscRawE_dy_sb->Scale(1./hscRawE_dy_sb->Integral());
  hs4_dy_sb->Scale(1./hs4_dy_sb->Integral());
  hr9_dy_sb->Scale(1./hr9_dy_sb->Integral());
  hsieip_dy_sb->Scale(1./hsieip_dy_sb->Integral());
  hsieie_dy_sb->Scale(1./hsieie_dy_sb->Integral());
  hesrr_dy_sb->Scale(1./hesrr_dy_sb->Integral());
  hEsOEn_dy_sb->Scale(1./hEsOEn_dy_sb->Integral());
  heta_dy_sb->Scale(1./heta_dy_sb->Integral());

  hetawidth_da_sb->Sumw2();
  hphiwidth_da_sb->Sumw2();
  hscRawE_da_sb->Sumw2();
  hs4_da_sb->Sumw2();
  hr9_da_sb->Sumw2();
  hsieip_da_sb->Sumw2();
  hsieie_da_sb->Sumw2();
  hesrr_da_sb->Sumw2();
  hEsOEn_da_sb->Sumw2();
  heta_da_sb->Sumw2();

  hetawidth_da_sb->Scale(1./hetawidth_da_sb->Integral());
  hphiwidth_da_sb->Scale(1./hphiwidth_da_sb->Integral());
  hscRawE_da_sb->Scale(1./hscRawE_da_sb->Integral());
  hs4_da_sb->Scale(1./hs4_da_sb->Integral());
  hr9_da_sb->Scale(1./hr9_da_sb->Integral());
  hsieip_da_sb->Scale(1./hsieip_da_sb->Integral());
  hsieie_da_sb->Scale(1./hsieie_da_sb->Integral());
  hesrr_da_sb->Scale(1./hesrr_da_sb->Integral());
  hEsOEn_da_sb->Scale(1./hEsOEn_da_sb->Integral());
  heta_da_sb->Scale(1./heta_da_sb->Integral());

  hetawidth_qcd_sb->Scale(1./hetawidth_qcd_sb->Integral());
  hphiwidth_qcd_sb->Scale(1./hphiwidth_qcd_sb->Integral());
  hscRawE_qcd_sb->Scale(1./hscRawE_qcd_sb->Integral());
  hs4_qcd_sb->Scale(1./hs4_qcd_sb->Integral());
  hr9_qcd_sb->Scale(1./hr9_qcd_sb->Integral());
  hsieip_qcd_sb->Scale(1./hsieip_qcd_sb->Integral());
  hsieie_qcd_sb->Scale(1./hsieie_qcd_sb->Integral());
  hesrr_qcd_sb->Scale(1./hesrr_qcd_sb->Integral());
  hEsOEn_qcd_sb->Scale(1./hEsOEn_qcd_sb->Integral());
  heta_qcd_sb->Scale(1./heta_qcd_sb->Integral());

  hetawidth_dy->SetLineColor(2);
  hphiwidth_dy->SetLineColor(2);
  hscRawE_dy->SetLineColor(2);
  hs4_dy->SetLineColor(2);
  hr9_dy->SetLineColor(2);
  hsieip_dy->SetLineColor(2);
  hsieie_dy->SetLineColor(2);
  hesrr_dy->SetLineColor(2);
  hEsOEn_dy->SetLineColor(2);
  heta_dy->SetLineColor(2);
  heta_dy->SetLineWidth(2);


  hetawidth_dy_sb->SetLineColor(4);
  hphiwidth_dy_sb->SetLineColor(4);
  hscRawE_dy_sb->SetLineColor(4);
  hs4_dy_sb->SetLineColor(4);
  hr9_dy_sb->SetLineColor(4);
  hsieip_dy_sb->SetLineColor(4);
  hsieie_dy_sb->SetLineColor(4);
  hesrr_dy_sb->SetLineColor(4);
  hEsOEn_dy_sb->SetLineColor(4);
  heta_dy_sb->SetLineColor(4);
  
  hetawidth_qcd_sb->SetLineColor(kGreen);
  hphiwidth_qcd_sb->SetLineColor(kGreen);
  hscRawE_qcd_sb->SetLineColor(kGreen);
  hs4_qcd_sb->SetLineColor(kGreen);
  hr9_qcd_sb->SetLineColor(kGreen);
  hsieip_qcd_sb->SetLineColor(kGreen);
  hsieie_qcd_sb->SetLineColor(kGreen);
  hesrr_qcd_sb->SetLineColor(kGreen);
  hEsOEn_qcd_sb->SetLineColor(kGreen);
  heta_qcd_sb->SetLineColor(kGreen);
  

  //set color for error histogram
  hetawidth_dy_err->SetLineColor(2);
  hphiwidth_dy_err->SetLineColor(2);
  hscRawE_dy_err->SetLineColor(2);
  hs4_dy_err->SetLineColor(2);
  hr9_dy_err->SetLineColor(2);
  hsieip_dy_err->SetLineColor(2);
  hsieie_dy_err->SetLineColor(2);
  hesrr_dy_err->SetLineColor(2);
  hEsOEn_dy_err->SetLineColor(2);
  heta_dy_err->SetLineColor(2);

  hetawidth_dy_sb_err->SetLineColor(4);
  hphiwidth_dy_sb_err->SetLineColor(4);
  hscRawE_dy_sb_err->SetLineColor(4);
  hs4_dy_sb_err->SetLineColor(4);
  hr9_dy_sb_err->SetLineColor(4);
  hsieip_dy_sb_err->SetLineColor(4);
  hsieie_dy_sb_err->SetLineColor(4);
  hesrr_dy_sb_err->SetLineColor(4);
  hEsOEn_dy_sb_err->SetLineColor(4);
  heta_dy_sb_err->SetLineColor(4);
  
  hetawidth_qcd_sb_err->SetLineColor(kGreen);
  hphiwidth_qcd_sb_err->SetLineColor(kGreen);
  hscRawE_qcd_sb_err->SetLineColor(kGreen);
  hs4_qcd_sb_err->SetLineColor(kGreen);
  hr9_qcd_sb_err->SetLineColor(kGreen);
  hsieip_qcd_sb_err->SetLineColor(kGreen);
  hsieie_qcd_sb_err->SetLineColor(kGreen);
  hesrr_qcd_sb_err->SetLineColor(kGreen);
  hEsOEn_qcd_sb_err->SetLineColor(kGreen);
  heta_qcd_sb_err->SetLineColor(kGreen);

  //set point is invisible
  hetawidth_dy_err->SetMarkerSize(0);
  hphiwidth_dy_err->SetMarkerSize(0);
  hscRawE_dy_err->SetMarkerSize(0);
  hs4_dy_err->SetMarkerSize(0);
  hr9_dy_err->SetMarkerSize(0);
  hsieip_dy_err->SetMarkerSize(0);
  hsieie_dy_err->SetMarkerSize(0);
  hesrr_dy_err->SetMarkerSize(0);
  hEsOEn_dy_err->SetMarkerSize(0);
  heta_dy_err->SetMarkerSize(0);

  hetawidth_dy_sb_err->SetMarkerSize(0);
  hphiwidth_dy_sb_err->SetMarkerSize(0);
  hscRawE_dy_sb_err->SetMarkerSize(0);
  hs4_dy_sb_err->SetMarkerSize(0);
  hr9_dy_sb_err->SetMarkerSize(0);
  hsieip_dy_sb_err->SetMarkerSize(0);
  hsieie_dy_sb_err->SetMarkerSize(0);
  hesrr_dy_sb_err->SetMarkerSize(0);
  hEsOEn_dy_sb_err->SetMarkerSize(0);
  heta_dy_sb_err->SetMarkerSize(0);
 
  hetawidth_qcd_sb_err->SetMarkerSize(0);
  hphiwidth_qcd_sb_err->SetMarkerSize(0);
  hscRawE_qcd_sb_err->SetMarkerSize(0);
  hs4_qcd_sb_err->SetMarkerSize(0);
  hr9_qcd_sb_err->SetMarkerSize(0);
  hsieip_qcd_sb_err->SetMarkerSize(0);
  hsieie_qcd_sb_err->SetMarkerSize(0);
  hesrr_qcd_sb_err->SetMarkerSize(0);
  hEsOEn_qcd_sb_err->SetMarkerSize(0);
  heta_qcd_sb_err->SetMarkerSize(0);
  

  c1 = new TCanvas("c1", "c1", 650, 550);
  c1->cd();
  c1->SetLogy();
  hetawidth_dy->GetYaxis()->SetTitle("a.u");
  hetawidth_dy->GetXaxis()->SetTitle("#eta^{width}_{SC}");
  hetawidth_dy->GetYaxis()->SetRangeUser(0.001, 1);
  hetawidth_dy->Draw("hist");
  hetawidth_dy_sb->Draw("histsame");
  hetawidth_qcd_sb->Draw("histsame");
  hetawidth_da_sb->Draw("same");
  hetawidth_dy_err->Draw("esames");
  hetawidth_dy_sb_err->Draw("esames");
  hetawidth_qcd_sb_err->Draw("esames");
  
  TLegend *lg1 = new TLegend(0.5, 0.65, 0.86, 0.88);
  lg1->AddEntry(hetawidth_dy, "signal region (DYJets)", "l");
  lg1->AddEntry(hetawidth_dy_sb, "side-band (DYJets)", "f");
  lg1->AddEntry(hetawidth_qcd_sb, "side-band (QCD)", "f");
  lg1->AddEntry(hetawidth_da_sb, "side-band (data)", "pe");
  lg1->SetTextFont(42);
  lg1->SetTextSize(0.04);
  lg1->SetBorderSize(0);
  lg1->Draw();
  TLatex tx;
  tx.SetNDC(kTRUE);
  tx.SetTextFont(62);
  tx.SetTextSize(0.04);
  tx.DrawLatex(0.6, 0.6, "EE 15 < p_{T}^{#gamma} < 20 GeV");

  c2 = new TCanvas("c2", "c2", 650, 550);
  c2->cd();
  c2->SetLogy();
  hphiwidth_dy->GetYaxis()->SetTitle("a.u");
  hphiwidth_dy->GetXaxis()->SetTitle("#phi^{width}_{SC}");
  hphiwidth_dy->GetYaxis()->SetRangeUser(0.001, 1);
  hphiwidth_dy->Draw("hist");
  hphiwidth_dy_sb->Draw("histsame");
  hphiwidth_qcd_sb->Draw("histsame");
  hphiwidth_da_sb->Draw("same");
  hphiwidth_dy_err->Draw("esames");
  hphiwidth_dy_sb_err->Draw("esames");
  hphiwidth_qcd_sb_err->Draw("esames");
  lg1->Draw();
  
  c3 = new TCanvas("c3", "c3", 650, 550);
  c3->cd();
  c3->SetLogy();
  hscRawE_dy->GetYaxis()->SetTitle("a.u");
  hscRawE_dy->GetXaxis()->SetTitle("E^{raw}_{SC} [GeV]");
  hscRawE_dy->GetYaxis()->SetRangeUser(0.001, 1);
  hscRawE_dy->Draw("hist");
  hscRawE_dy_sb->Draw("histsame");
  hscRawE_qcd_sb->Draw("histsame");
  hscRawE_da_sb->Draw("same");
  hscRawE_dy_err->Draw("esames");
  hscRawE_dy_sb_err->Draw("esames");
  hscRawE_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c4 = new TCanvas("c4", "c4", 650, 550);
  c4->cd();
  c4->SetLogy();
  hs4_dy->GetYaxis()->SetTitle("a.u");
  hs4_dy->GetXaxis()->SetTitle("S4");
  hs4_dy->GetYaxis()->SetRangeUser(0.001, 6);
  hs4_dy->Draw("hist");
  hs4_dy_sb->Draw("histsame");
  hs4_qcd_sb->Draw("histsame");
  hs4_da_sb->Draw("same");
  hs4_dy_err->Draw("esames");
  hs4_dy_sb_err->Draw("esames");
  hs4_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c5 = new TCanvas("c5", "c5", 650, 550);
  c5->cd();
  c5->SetLogy();
  hr9_dy->GetYaxis()->SetTitle("a.u");
  hr9_dy->GetXaxis()->SetTitle("R9");
  hr9_dy->GetYaxis()->SetRangeUser(0.001, 6);
  hr9_dy->Draw("hist");
  hr9_dy_sb->Draw("histsame");
  hr9_qcd_sb->Draw("histsame");
  hr9_da_sb->Draw("same");
  hr9_dy_err->Draw("esames");
  hr9_dy_sb_err->Draw("esames");
  hr9_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c6 = new TCanvas("c6", "c6", 650, 550);
  c6->cd();
  c6->SetLogy();
  hsieip_dy->GetYaxis()->SetTitle("a.u");
  hsieip_dy->GetXaxis()->SetTitle("#sigma_{i#etai#phi}");
  hsieip_dy->GetYaxis()->SetRangeUser(0.001, 4);
  hsieip_dy->Draw("hist");
  hsieip_dy_sb->Draw("histsame");
  hsieip_qcd_sb->Draw("histsame");
  hsieip_da_sb->Draw("same");
  hsieip_dy_err->Draw("esames");
  hsieip_dy_sb_err->Draw("esames");
  hsieip_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c7 = new TCanvas("c7", "c7", 650, 550);
  c7->cd();
  c7->SetLogy();
  hsieie_dy->GetYaxis()->SetTitle("a.u");
  hsieie_dy->GetXaxis()->SetTitle("#sigma_{i#etai#eta}");
  hsieie_dy->GetYaxis()->SetRangeUser(0.001, 4);
  hsieie_dy->Draw("hist");
  hsieie_dy_sb->Draw("histsame");
  hsieie_qcd_sb->Draw("histsame");
  hsieie_da_sb->Draw("same");
  hsieie_dy_err->Draw("esames");
  hsieie_dy_sb_err->Draw("esames");
  hsieie_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c8 = new TCanvas("c8", "c8", 650, 550);
  c8->cd();
  c8->SetLogy();
  hesrr_dy->GetYaxis()->SetTitle("a.u");
  hesrr_dy->GetXaxis()->SetTitle("#sigma^{ES}");
  hesrr_dy->GetYaxis()->SetRangeUser(0.001, 1);
  hesrr_dy->Draw("hist");
  hesrr_dy_sb->Draw("histsame");
  hesrr_qcd_sb->Draw("histsame");
  hesrr_da_sb->Draw("same");
  hesrr_dy_err->Draw("esames");
  hesrr_dy_sb_err->Draw("esames");
  hesrr_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c9 = new TCanvas("c9", "c9", 650, 550);
  c9->cd();
  c9->SetLogy();
  hEsOEn_dy->GetYaxis()->SetTitle("a.u");
  hEsOEn_dy->GetXaxis()->SetTitle("ES energy/ Raw energy");
  hEsOEn_dy->GetYaxis()->SetRangeUser(0.001, 1);
  hEsOEn_dy->Draw("hist");
  hEsOEn_dy_sb->Draw("histsame");
  hEsOEn_qcd_sb->Draw("histsame");
  hEsOEn_da_sb->Draw("same");
  hEsOEn_dy_err->Draw("esames");
  hEsOEn_dy_sb_err->Draw("esames");
  hEsOEn_qcd_sb_err->Draw("esames");
  lg1->Draw();

  c10 = new TCanvas("c10", "c10", 650, 550);
  c10->cd();
  c10->SetLogy();
  heta_dy->GetYaxis()->SetTitle("a.u");
  heta_dy->GetXaxis()->SetTitle("#eta");
  heta_dy->GetYaxis()->SetRangeUser(0.001, 1);
  heta_dy->Draw("hist");
  heta_dy_sb->Draw("histsame");
  //heta_qcd_sb->Draw("histsame");
  heta_da_sb->Draw("same");
  heta_dy_err->Draw("esames");
  heta_dy_sb_err->Draw("esames");
  //heta_qcd_sb_err->Draw("esames");
  TLegend *lg = new TLegend(0.35, 0.28, 0.68, 0.56);
  lg->AddEntry(heta_dy, "signal region (DYJets)", "l");
  lg->AddEntry(heta_dy_sb, "side-band (DYJets)", "f");
  //lg->AddEntry(heta_qcd_sb, "side-band (QCD)", "f");
  lg->AddEntry(heta_da_sb, "side-band (data)", "pe");
  lg->SetTextFont(42);
  lg->SetTextSize(0.04);
  lg->SetBorderSize(0);
  lg->Draw();

  TH1F *hratio = (TH1F*) heta_dy->Clone();
  hratio->Divide(heta_dy_sb);
  hratio->Print("all");

  cn = new TCanvas("cn", "cn", 650, 550);
  cn->cd();
  hratio->GetYaxis()->SetTitle("SR/SB");
  hratio->GetXaxis()->SetTitle("#eta");
  hratio->GetYaxis()->SetRangeUser(0, 3);
  hratio->SetLineColor(1);
  hratio->Draw();
  hratio->Fit("expo", "", "", 1.5, 2.5);

  c1->SaveAs("plots/etawidth_SR_SB_Pt15to20_ele_EE_mc_data_etawei.pdf");
  //c2->SaveAs("plots/phiwidth_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  c3->SaveAs("plots/rawEn_SR_SB_Pt15to20_ele_EE_mc_data_etawei.pdf");
  //c4->SaveAs("plots/s4_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c5->SaveAs("plots/r9_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c6->SaveAs("plots/sieip_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c7->SaveAs("plots/sieie_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c8->SaveAs("plots/ESsigma_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c9->SaveAs("plots/ESOEn_SR_SB_Pt15to20_ele_EE_mc_data.pdf");
  //c10->SaveAs("plots/eta_SR_SB_Pt15to20_ele_EE_mc_data_etawei.pdf");



}
