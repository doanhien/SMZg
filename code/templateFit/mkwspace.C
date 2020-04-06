void mkwspace(bool ele = true, bool barrel = true, TString cat = "pt", float minpt = 15, float maxpt = 20){
  // As usual, load the combine library to get access to the RooParametricHist
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  

  //TString dir = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/";
  TString dir = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/exclusive/";
  //TString dir = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/signalTemplate_2sigma/";
  //TString dir = "/afs/cern.ch/work/t/tdoan/analysis/SMZg/templatefit/PhoPt20/SB_EB7to13_EE6to14/";
  TString fname = dir;

  if (cat == "pt") {
    //fname += Form("Ptg/datacard/data_Pt%dto%d", (int) minpt, (int) maxpt);
    fname += Form("Ptllg/datacard/data_Ptllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) fname += "_EB";
    else fname += "_EE";
    if (ele) fname += "_ele.root";
    else fname += "_mu.root";
  }
  else {
    fname += Form("Mllg/datacard/data_Mllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) fname += "_EB";
    else fname += "_EE";

    if (ele) fname += "_ele.root";
    else fname += "_mu.root";
  }


  TFile *fin = new TFile (fname, "read");
  TH1F *hdata = (TH1F*)fin->Get("data_obs");
  TH1F *hsig  = (TH1F*)fin->Get("signal");
  TH1F *hsig_sigmaUp = (TH1F*)fin->Get("signal_sigmaUp");
  TH1F *hsig_sigmaDown = (TH1F*)fin->Get("signal_sigmaDown");

  TH1F *hbkg = (TH1F*)fin->Get("background");
  

  TString fitname = dir;

  if (cat == "pt") {
    //fitname += Form("Ptg/1stFit/fitDiagnostics_Pt%dto%d", (int) minpt, (int) maxpt);
    fitname += Form("Ptllg/1stFit/fitDiagnostics_Ptllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) fitname += "_EB";
    else fitname += "_EE";
    if (ele) fitname += "_ele.root";
    else fitname += "_mu.root";
  }
  else {
    fitname += Form("Mllg/1stFit/fitDiagnostics_Mllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) fitname += "_EB";
    else fitname += "_EE";
    if (ele) fitname += "_ele.root";
    else fitname += "_mu.root";
  }

  TFile *ffit = new TFile(fitname, "read");

  RooArgSet *norm_fit_s = (RooArgSet*) ffit->Get("norm_fit_s");
  RooRealVar* parfit_sig = (RooRealVar*) norm_fit_s->find("Zg/signal");
  RooRealVar* parfit_bkg = (RooRealVar*) norm_fit_s->find("Zg/background");
  float  signalYield = parfit_sig->getValV();
  float  bkgYield = parfit_bkg->getValV();

  cout << "signal yield from rateParam fit: " << signalYield << endl;

  hsig->Scale(signalYield/hsig->Integral());
  hsig_sigmaUp->Scale(signalYield/hsig_sigmaUp->Integral());
  hsig_sigmaDown->Scale(signalYield/hsig_sigmaDown->Integral());


  hdata->SetName("data_obs");
  hsig->SetName("signal");
  hsig_sigmaUp->SetName("signal_sigmaUp");
  hsig_sigmaDown->SetName("signal_sigmaDown");
  hbkg->SetName("background_ref");

  float bkg_template_err[10];
  for(int ii=0; ii<10; ii++){
    bkg_template_err[ii] = TMath::Sqrt(hbkg->GetBinContent(ii+1))/hbkg->GetBinContent(ii+1);
  }


  TH1F *hsideband = (TH1F*)hbkg->Clone();
  if(bkgYield>=0 && bkgYield<1.) bkgYield=1.;

  hsideband->Scale(bkgYield/hsideband->Integral());  
  hsideband->SetName("background");

  //hsig->Scale(1./hsig->Integral());
  //hsideband->Scale(1./hsideband->Integral());

  printf("observation %d \n", int(hdata->Integral()));
  printf("rate            %.1f      %.1f \n", hsig->Integral(), hsideband->Integral());


  // Output file and workspace 
  TString outdir = dir;

  //if (cat == "pt") outdir += "Ptg/2ndFit/datacard/";
  if (cat == "pt") outdir += "Ptllg/2ndFit/datacard/";
  else outdir += "Mllg/2ndFit/datacard/";

  const int dir_err = system("mkdir -p " + outdir);
  if (-1 == dir_err)
    {
      printf("Error creating directory!n");
      exit(1);
    }


  TString outname = outdir;
  if (cat == "pt") {
    //outname += Form("Pt%dto%d", (int) minpt, (int) maxpt);
    outname += Form("Ptllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) outname += "_EB";
    else outname += "_EE";
    if (ele) outname += "_ele_ws.root";
    else outname += "_mu_ws.root";
  }
  else {
    outname += Form("Mllg%dto%d", (int) minpt, (int) maxpt);
    if (barrel) outname += "_EB";
    else outname += "_EE";
    if (ele) outname += "_ele_ws.root";
    else outname += "_mu_ws.root";
  }

  TFile *fOut = new TFile(outname,"RECREATE");
  RooWorkspace wspace("wspace","wspace");

  // A fit in BDT output
  RooRealVar bdt("BDT","BDT",-1., 1.);
  RooArgList vars(bdt);

  RooDataHist data("data_obs","data_obs", vars, hdata);
  RooDataHist signal("signal","signal", vars, hsig);
  RooDataHist signal_sigmaUp("signal_sigmaUp","signal_sigmaUp", vars, hsig_sigmaUp);
  RooDataHist signal_sigmaDown("signal_sigmaDown","signal_sigmaDown", vars, hsig_sigmaDown);

  RooRealVar bin01("bkg_bin01","bkg bin 01", hsideband->GetBinContent(1), 0., hsideband->Integral());
  RooRealVar bin02("bkg_bin02","bkg bin 02", hsideband->GetBinContent(2), 0., hsideband->Integral());
  RooRealVar bin03("bkg_bin03","bkg bin 03", hsideband->GetBinContent(3), 0., hsideband->Integral());
  RooRealVar bin04("bkg_bin04","bkg bin 04", hsideband->GetBinContent(4), 0., hsideband->Integral());
  RooRealVar bin05("bkg_bin05","bkg bin 05", hsideband->GetBinContent(5), 0., hsideband->Integral());
  RooRealVar bin06("bkg_bin06","bkg bin 06", hsideband->GetBinContent(6), 0., hsideband->Integral());
  RooRealVar bin07("bkg_bin07","bkg bin 07", hsideband->GetBinContent(7), 0., hsideband->Integral());
  RooRealVar bin08("bkg_bin08","bkg bin 08", hsideband->GetBinContent(8), 0., hsideband->Integral());
  RooRealVar bin09("bkg_bin09","bkg bin 09", hsideband->GetBinContent(9), 0., hsideband->Integral());
  RooRealVar bin10("bkg_bin10","bkg bin 10", hsideband->GetBinContent(10), 0., hsideband->Integral());

  for(int ii=0; ii<10; ii++){
    printf("bkg_bin%02d  param  %.1f  %.1f\n", ii, hsideband->GetBinContent(ii+1), hsideband->GetBinContent(ii+1)*bkg_template_err[ii]);
  }

  RooArgList bkg_bins;
  bkg_bins.add(bin01);
  bkg_bins.add(bin02);
  bkg_bins.add(bin03);
  bkg_bins.add(bin04);
  bkg_bins.add(bin05);
  bkg_bins.add(bin06);
  bkg_bins.add(bin07);
  bkg_bins.add(bin08);
  bkg_bins.add(bin09);
  bkg_bins.add(bin10);
  
  RooParametricHist bkg("background","background", bdt, bkg_bins, *hsideband);
  RooAddition bkg_norm("background_norm","Total Number of events from background in signal region", bkg_bins);

  //parameters for signal
  RooRealVar bin01_sig("sig_bin01","sig bin 01", hsig->GetBinContent(1), 0., hsig->Integral());
  RooRealVar bin02_sig("sig_bin02","sig bin 02", hsig->GetBinContent(2), 0., hsig->Integral());
  RooRealVar bin03_sig("sig_bin03","sig bin 03", hsig->GetBinContent(3), 0., hsig->Integral());
  RooRealVar bin04_sig("sig_bin04","sig bin 04", hsig->GetBinContent(4), 0., hsig->Integral());
  RooRealVar bin05_sig("sig_bin05","sig bin 05", hsig->GetBinContent(5), 0., hsig->Integral());
  RooRealVar bin06_sig("sig_bin06","sig bin 06", hsig->GetBinContent(6), 0., hsig->Integral());
  RooRealVar bin07_sig("sig_bin07","sig bin 07", hsig->GetBinContent(7), 0., hsig->Integral());
  RooRealVar bin08_sig("sig_bin08","sig bin 08", hsig->GetBinContent(8), 0., hsig->Integral());
  RooRealVar bin09_sig("sig_bin09","sig bin 09", hsig->GetBinContent(9), 0., hsig->Integral());
  RooRealVar bin10_sig("sig_bin10","sig bin 10", hsig->GetBinContent(10), 0., hsig->Integral());

  RooArgList sig_bins;
  sig_bins.add(bin01_sig);
  sig_bins.add(bin02_sig);
  sig_bins.add(bin03_sig);
  sig_bins.add(bin04_sig);
  sig_bins.add(bin05_sig);
  sig_bins.add(bin06_sig);
  sig_bins.add(bin07_sig);
  sig_bins.add(bin08_sig);
  sig_bins.add(bin09_sig);
  sig_bins.add(bin10_sig);
  
  RooParametricHist sig("signal","signal", bdt, sig_bins, *hsig);
  RooAddition sig_norm("signal_norm","Total Number of signal", sig_bins);


 // import the pdfs
  wspace.import(data);
  wspace.import(signal);
  wspace.import(signal_sigmaUp);
  wspace.import(signal_sigmaDown);

  wspace.import(bkg);
  wspace.import(bkg_norm,RooFit::RecycleConflictNodes());
  //wspace.import(sig);
  //wspace.import(sig_norm,RooFit::RecycleConflictNodes());


  wspace.Write();

  hdata->Write();
  hsig->Write();
  hsig_sigmaUp->Write();
  hsig_sigmaDown->Write();
  hbkg->Write();
  hsideband->Write();

  //this for datacard 
  TString sufix;
  if (cat == "pt" )
    //sufix= Form("_Pt%dto%d", (int) minpt, (int) maxpt);
    sufix= Form("_Ptllg%dto%d", (int) minpt, (int) maxpt);
  else
    sufix= Form("_Mllg%dto%d", (int) minpt, (int) maxpt);

  //  if ( cat == "pt" ) {
    if (barrel) sufix += "_EB";
    else sufix += "_EE";
    //}
  if (ele) sufix += "_ele";
  else sufix += "_mu";

  TString fnamecard = outdir;
  fnamecard += "datacard";
  fnamecard += sufix;
  fnamecard += "_ws.txt";

  ofstream fcard;
  fcard.open(fnamecard);
  cout << "outname name: " << fnamecard.Data() << endl;

  cout << "write data card" << endl;

  fcard << "#### range:" << (int) minpt <<"-" << (int) maxpt << " GeV";
  if (cat == "pt") {
    if (barrel) fcard << " in EB,";
    else fcard << "EE,";
  }
  if (ele) fcard << "ele channel" << endl;
  else fcard << "muon channel" << endl;

  fcard << "------------------" << endl;
  fcard << "imax *" << endl;
  fcard << "jmax *" << endl;
  fcard << "kmax *" << endl;
  fcard << "------------------" << endl;

  fcard << "shapes * * " << outname.Data() << " wspace:$PROCESS" << " wspace:$PROCESS_$SYSTEMATIC" << endl;
  fcard << "------------------" << endl;
  fcard << "bin Zg" << endl;
  fcard << "observation " << -1 << endl;
  fcard << "----------------------------------------------" << endl;
  fcard << "bin \t \t Zg \t \t Zg" << endl;
  fcard << "process \t signal \t background" << endl;
  fcard << "process \t 0 \t \t 1" << endl;
  fcard << "rate \t \t -1 \t  1" << endl;
  fcard << "----------------------------------------------" << endl;
  fcard << "sigma \t shape \t 1.0 \t-" << endl;

  TH1 *hbkg_rateParaFit = (TH1D*) ffit->Get("shapes_fit_s/total_background");

  for (int i = 1; i <= hsideband->GetNbinsX(); i++) {
    float tmpErr = hsideband->GetBinError(i);
    if(tmpErr < 0.1) {
      if (hsideband->GetBinContent(i)>=0.1)  tmpErr = 0.1;
      else tmpErr = hsideband->GetBinContent(i) ;
    }
    //if(tmpErr < 0.02) tmpErr = 0.02;

    float bkgyield = hsideband->GetBinContent(i);
    if ( bkgyield < 1e-6) {
      bkgyield = 1e-6;
      tmpErr = 1e-6;
    }

    if (i < 10) fcard << "bkg_bin0" << i << "  param  " << bkgyield << "   " << tmpErr << endl;
    else fcard << "bkg_bin" << i << "  param  "  << bkgyield << "   " << tmpErr << endl;
  }


  fcard.close();

  cout << "done!!!!!!" << endl;


  // f->Close();  
  fin->Close();
  fOut->Close();
  

}
