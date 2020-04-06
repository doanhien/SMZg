{
  gSystem->SetBuildDir("tmpdir", kTRUE);
  gROOT->LoadMacro("templateMaker.C+");   

  Float_t XS[4] = {1, 117.864, 5943.2, 1};
  TString sample[9] = {
    //"job_spring16_gjet_pt15to6000.root",
    "job_summer16_WGToLNuG_amcatnlo_ext3.root",
    "job_summer16_GJets_Pt_15To6000_TuneCUETP8M1_pythia8.root",
    "job_summer16_WJetsToLNu_amcatnlo_ext2.root",
    "job_summer16_QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8.root",
    "job_summer16_QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8.root",
    "job_summer16_QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8.root",
    "job_summer16_QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8.root",
    "job_summer16_QCD_Pt_80to120_TuneCUETP8M1_pythia8_ext2.root",
    "JetHT_Run2016B_ReMiniAOD_03Feb.root",
    //"job_summer16_GJets_Pt_15To6000_TuneCUETP8M1_pythia8.root",
    //"job_summer16_Zg_aMCatNLO.root",
    //"job_summer16_DYJetsToLL_m50_aMCatNLO.root"
    //"job_spring16_Wg_MG.root",
    //"job_spring16_gjet_pt20_MGG_40to80.root",
    //"job_spring16_gjet_pt20to40_MGG_80toInf.root",
    //"job_spring16_gjet_pt40_MGG_80toInf.root"
  };
  Float_t Lumi           = 35900; // pb^-1
  Int_t aMCatNLO[9]      = {1, 0, 1, 0, 0, 0, 0, 0, 0};

  for (Int_t i=3; i<6; ++i) {

    //TString inpath  = "/data6/ggNtuples/V08_00_11_01/"+sample[i];
    //TString inpath  = "/data6/ggNtuples/V08_00_26_05/skim/"+sample[i];
    TString inpath;
    inpath  = "/data1/ggNtuples/V826/mc/skim/" +sample[i];
    if (i==8) inpath = "/data1/ggNtuples/V826/data/" + sample[i];
    
    const char* inpaths[]  = {inpath.Data()};
    //TString outname        = "bkgTemplate_summer16_420_LooseSieie";
    TString outname        = "bkgTemplate_summer16_420_BDT_Up6000_5VarCorr";
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
