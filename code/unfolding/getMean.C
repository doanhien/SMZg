#include "TFile.h"

void getMean(TString filename = "EGM_ER_ES_Unc/ele_resol_unc_ptllg_ele_EB.root") {

  cout << "open file :" << filename << endl;
  TFile *f = new TFile(filename, "read");
  TTree *t1 = (TTree*) f->Get("outtree");

  cout << setprecision(3) << std::fixed << endl;

  int nbins = t1->GetMaximum("nbins");
  cout << "nbins: " << nbins << endl;

  cout << "\n -------- uncertainty" << endl;
  for (int i = 0; i < nbins; i++) {
    TH1F *h1 = new TH1F("h1", "h1", 200, 0, 200);
    t1->Draw(Form("yield[%d] >> h1", i), "", "goff");
    cout << h1->GetMean() << " ";
    h1->Delete();
  }

  cout << "\n done" << endl;
}

