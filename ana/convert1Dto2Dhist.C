#include "TFile.h"
#include "TH1.h"
#include "TH2.h"


void convert1Dto2Dhist() {

  TFile *f1 = new TFile("external/SFs_pt_Mu8Leg_Eta0to09_PtBelow200.root", "READ");

  TH1F *h1 = (TH1F*) f1->Get("scale_factor");

  TFile *f[4];
  TH1F *h[4];

  f[0] = new TFile("external/SFs_pt_Mu8Leg_Eta0to09_PtBelow200.root", "READ");
  f[1] = new TFile("external/SFs_pt_Mu8Leg_Eta09to12_PtBelow200.root", "READ");
  f[2] = new TFile("external/SFs_pt_Mu8Leg_Eta12to21_PtBelow200.root", "READ");
  f[3] = new TFile("external/SFs_pt_Mu8Leg_Eta21to24_PtBelow200.root", "READ");

  for (int ifi = 0; ifi < 4; ifi++) {
    h[ifi] = (TH1F*) f[ifi]->Get("scale_factor");
  }


  int NPtbin = h1->GetNbinsX();
  cout << "number of pt bin: " << NPtbin << endl;
  float Pt[NPtbin+1];

  for (int i = 0; i < NPtbin; i ++) {
    Pt[i] = h1->GetXaxis()->GetBinLowEdge(i+1);
    //cout << "Pt: " << Pt[i] << endl;
  }
  Pt[NPtbin] = h1->GetXaxis()->GetBinUpEdge(NPtbin);
  //cout << Pt[NPtbin+1] << endl;

  for (int i =0; i <= NPtbin; i++) {
    cout << i << "\tPt: " << Pt[i] << endl; 
  }

  const int NEtabin = 4;
  Float_t etabin[NEtabin+1] = {0.0, 0.9, 1.2, 2.1, 2.4};
  TH2D *scale_factor = new TH2D("scale_factor", "scale_factor", NEtabin, etabin, NPtbin, Pt);

  for (int ix = 1; ix <= NEtabin; ix++) {
    for (int iy = 1; iy <= NPtbin; iy++) {
      float sf = h[ix-1]->GetBinContent(iy, ix);
      float errsf = h[ix-1]->GetBinError(iy, ix);
      cout << "SF: " << sf << endl;
      scale_factor->SetBinContent(ix, iy, sf);
      scale_factor->SetBinError(ix, iy, errsf);
    }
  }

  //c1 = new TCanvas("c1", "c1", 650, 550);
  //c1->cd();
  //scale_factor->Draw("colz text");

  TFile *fout = new TFile("external/SFs_Leg2_Mu8_PtBelow200.root", "recreate");
  fout->cd();
  scale_factor->Write();
  fout->Write();
  fout->Close();




}
