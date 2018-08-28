{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->LoadMacro("xAna.C+");

  int iWP = 0;

  //runData(0, iWP);
  runMC(1, iWP);

}

void runData(bool mc = 0, int iWP = 0) {

  TString sample[9] = {
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016B_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016C_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016D_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016E_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016F_FebReminiAOD1.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016F_FebReminiAOD2.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016G_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016H_FebReminiAODv2.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016H_FebReminiAODv3.root"

  };

  /*
  //muon dataset
  TString sample[9] = {
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016B_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016C_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016D_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016E_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016F_FebReminiAOD1.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016F_FebReminiAOD2.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016G_FebReminiAOD.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016H_FebReminiAODv2.root",
    "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleMu_Run2016H_FebReminiAODv3.root",
  };
  */

  for (Int_t i=0; i<9; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    xAna(inpaths, "out.root",  mc, iWP);
  }
}

void runMC(bool mc = 1, int iWP = 0) {


  TString sample[10] = {
    "/data6/ggNtuples/V08_00_26_05/job_summer16_ZGToLLG_5f_lhe_FixedEta_FixedPhi_beforefsr_skim.root",
    //"/data6/ggNtuples/V08_00_26_05/skim/job_summer16_Zg_aMCatNLO.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_DYJetsToLL_m50_aMCatNLO.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_TTTo2L2Nu_powheg.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WWTo2L2Nu.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WWToLNuQQ*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WZTo3LNu.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_WZTo2L2Q.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo2L2Nu*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo2L2Q.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_ZZTo4L.root",
  };


  for (Int_t i=0; i<1; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    xAna(inpaths, "out.root", mc, iWP);
  }

}
