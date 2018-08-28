{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  //gROOT->LoadMacro("SigTemplate.C+");   
  gROOT->LoadMacro("templateMaker.C+");   

  Float_t XS[12] = {1, 117.864, 5943.2, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  TString sample[12] = {
    "job_spring16_gjet_pt15to6000.root",
    "job_summer16_Zg_aMCatNLO.root",
    "job_summer16_DYJetsToLL_m50_aMCatNLO.root",
    "job_DoubleEG_Run2016B_FebReminiAOD.root",
    "job_DoubleEG_Run2016C_FebReminiAOD.root",
    "job_DoubleEG_Run2016D_FebReminiAOD.root",
    "job_DoubleEG_Run2016E_FebReminiAOD.root",
    "job_DoubleEG_Run2016F_FebReminiAOD1.root", 
    "job_DoubleEG_Run2016F_FebReminiAOD2.root",
    "job_DoubleEG_Run2016G_FebReminiAOD.root",
    "job_DoubleEG_Run2016H_FebReminiAODv2.root",
    "job_DoubleEG_Run2016H_FebReminiAODv3.root"
    
    /*
    "job_DoubleMu_Run2016B_FebReminiAOD.root",
    "job_DoubleMu_Run2016C_FebReminiAOD.root",
    "job_DoubleMu_Run2016D_FebReminiAOD.root",
    "job_DoubleMu_Run2016E_FebReminiAOD.root",
    "job_DoubleMu_Run2016F_FebReminiAOD1.root", 
    "job_DoubleMu_Run2016F_FebReminiAOD2.root",
    "job_DoubleMu_Run2016G_FebReminiAOD.root",
    "job_DoubleMu_Run2016H_FebReminiAODv2.root",
    "job_DoubleMu_Run2016H_FebReminiAODv3.root"
    */
    //"job_spring16_Wg_MG.root",
    //"job_spring16_gjet_pt20_MGG_40to80.root",
    //"job_spring16_gjet_pt20to40_MGG_80toInf.root",
    //"job_spring16_gjet_pt40_MGG_80toInf.root"
  };
  Float_t Lumi           = 35900; // pb^-1
  Int_t aMCatNLO[12]      = {0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  for (Int_t i=1; i<3; ++i) {

    TString inpath;
    if (i == 0) inpath = "/data6/ggNtuples/V08_00_11_01/";
    else if ( i==1 || i ==2) inpath  = "/data6/ggNtuples/V08_00_26_05/skim/";
    else inpath  = "/data1/ggNtuples/V08_00_26_05/skim/";
    inpath += sample[i];

    //TString inpath  = "/data6/ggNtuples/V08_00_11_01/"+sample[i];
    //TString inpath  = "/data6/ggNtuples/V08_00_26_05/skim/"+sample[i];
    //TString inpath  = "/data1/ggNtuples/V08_00_26_05/skim/"+sample[i];
    
    const char* inpaths[]  = {inpath.Data()};
    TString outname        = "bkgTemplate_summer16_420_LooseSieie";
    outname = "output/" + outname + "_" +sample[i];
    cout<<outname<<endl;
    const char* outpath    = outname.Data();
    Float_t xsection       = XS[i];
    Int_t   negativeWeight = aMCatNLO[i];

    cout<<"processing : "<<inpath<<endl;
    //inpaths, npaths, outpath, xs=1, lumi=1, aMCatNLO=0
    templateMaker(inpaths, 1, outpath, xsection, Lumi, negativeWeight);

  }

}
