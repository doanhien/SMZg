#include "Riostream.h"

void draw() {

  ifstream in_1, in_2, in_3, in_4, in_5;

  in_1.open("");
  in_2.open("");
  in_3.open("");
  in_4.open("");
  in_5.open("");

  Float_t y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12;
  Int_t nlines = 0;
  TFile *f = new TFile("yield_Rochester_stat_SB_EB7to13_EE6to14_Bkg_Subtract_PhoPt20_excl.root","RECREATE");
  TH1F *h1 = new TH1F("h1","y distribution",100,-4,4);
  TNtuple *ntuple_barrel = new TNtuple("ntuple_barrel","data from ascii file","y1:y2:y3:y4:y5:y6:y7:y8:y9:y10:y11"); //for pt
  TNtuple *ntuple_endcap = new TNtuple("ntuple_endcap","data from ascii file","y1:y2:y3:y4:y5:y6:y7:y8:y9:y10:y11"); //for pt
  TNtuple *ntuple_mllg_barrel = new TNtuple("ntuple_mllg_barrel","data from ascii file","y1:y2:y3:y4:y5:y6:y7:y8:y9:y10:y11");     // for mllg
  TNtuple *ntuple_mllg_endcap = new TNtuple("ntuple_mllg_endcap","data from ascii file","y1:y2:y3:y4:y5:y6:y7:y8:y9:y10:y11");     // for mllg
  TNtuple *ntuple_njet = new TNtuple("ntuple_njet","data from ascii file","y1:y2:y3:y4"); //for njet

  while (1) {
    cout << "fill ntuple" << endl;
    in_1 >> y1 >> y2 >> y3 >> y4 >> y5 >> y6 >> y7 >> y8 >> y9 >> y10 >> y11; //for pt
    if (!in_1.good()) break;
    ntuple_barrel->Fill(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11);
    cout << "fill endcap " << endl;
    in_2 >> y1 >> y2 >> y3 >> y4 >> y5 >> y6 >> y7 >> y8 >> y9 >> y10 >> y11; //for pt
    if (!in_2.good()) break;
    ntuple_endcap->Fill(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11);

    in_3 >> y1 >> y2 >> y3 >> y4 >> y5 >> y6 >> y7 >> y8 >> y9 >> y10 >> y11; //for mllg
    if (!in_3.good()) break;
    ntuple_mllg_barrel->Fill(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11);

    cout << "fill EE for mllg" << endl;
    in_4 >> y1 >> y2 >> y3 >> y4 >> y5 >> y6 >> y7 >> y8 >> y9 >> y10 >> y11; //for mllg 
    if (!in_4.good()) break;
    ntuple_mllg_endcap->Fill(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11);

    cout << "filling njet" << endl;
    in_5 >> y1 >> y2 >> y3 >> y4; //for njet
    if (!in_5.good()) break;
    cout << "done njet" << endl;

    if (nlines < 5) printf("y1=%8f, y2=%8f, y3=%8f\n",y1,y2,y3);
    h1->Fill(y1);
    ntuple_njet->Fill(y1, y2, y3, y4);
    cout << "done filling" << endl;

    nlines++;
  }
  printf(" found %d points\n",nlines);

  in_1.close();
  in_2.close();
  in_3.close();
  in_4.close();
  in_5.close();

  f->Write();
}


}
