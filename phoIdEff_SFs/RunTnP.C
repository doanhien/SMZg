//void RunTnP(bool data = 0, int nfile = 1)
{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->LoadMacro("anaTnP_v2.C+");

  //runData(0);
  runMC(1);

}

void runData(bool mc = 0) {
  //void runData( int nfile = 0, bool mc = 0) {

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

    //TString inpaths;
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016B_FebReminiAOD.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016C_FebReminiAOD.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016D_FebReminiAOD.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016E_FebReminiAOD.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016F_FebReminiAOD1.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016F_FebReminiAOD2.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016G_FebReminiAOD.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016H_FebReminiAODv2.root";
  //inpaths = "/data1/ggNtuples/V08_00_26_05/skim/job_DoubleEG_Run2016H_FebReminiAODv3.root";


  for (Int_t i=0; i<9; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    anaTnP_v2(inpaths, "out.root", mc);
  }

}

//void runMC(bool mc = 1, bool ele = 1, int iWP = 0) {
void runMC(bool mc = 1) {

  TString inpaths;
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_XZg_W0p014_m1500/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_XZg_W1p4_m500/ggtree_mc*.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_Zg_aMCatNLO.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_Zg_pt130.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_DYJetsToLL_m50_aMCatNLO.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_WZTo3LNu.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_ZZTo4L.root";
  //inpaths = "/data6/ggNtuples/V08_00_11_01/job_spring16_TT_powheg_ext3.root";

  inpaths = "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_DYJetsToLL_m50_aMCatNLO.root";
  //inpaths = "/data6/ggNtuples/V08_00_24_00/summer16_DYJetsToLL_m50_MG/ggtree_mc*.root";

  cout << "Process sample: " << inpaths << endl;
  anaTnP_v2(inpaths, "out.root", mc);

}
