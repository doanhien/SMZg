#include "anaTnP_v2.C"

void runMC(bool);
void runData(bool);

void RunTnP(bool data = 0, int nfile = 1)
{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  //gROOT->LoadMacro("anaTnP_v2.C+");

  //runData(0);
  runMC(1);

}

void runData(bool mc = 0) {
  //void runData( int nfile = 0, bool mc = 0) {

  TString sample[9] = {
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016B_FebReminiAOD/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016C_FebReminiAOD/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016D_FebReminiAOD/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016E_FebReminiAOD/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016F_FebReminiAOD1/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016F_FebReminiAOD2/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016G_FebReminiAOD/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016H_FebReminiAOD2/ggtree_data*.root",
    "/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016H_FebReminiAOD3/ggtree_data*.root",

  };

  for (Int_t i=0; i<9; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    anaTnP_v2(inpaths, "out.root", mc);
  }

}


void runMC(bool mc = 1) {

  TString sample[2] = {
    "/data6/ggNtuples/V09_04_13_04/job_summer16_DYJetsToLL_m50_aMCatNLO_ext2/ggtree*.root",
    "/data6/ggNtuples/V08_00_26_05/skim/job_summer16_TTTo2L2Nu_powheg.root"
  };

  for (Int_t i=0; i<1; ++i) {
    TString inpaths = sample[i].Data();
    cout << "Process sample: " << inpaths << endl;
    anaTnP_v2(inpaths, "out.root", mc);
  }

}
