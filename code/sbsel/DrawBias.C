#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void DrawBias (bool ele = true, string cat = "EB", float minpt = 15., float maxpt = 20.) {


  gStyle->SetOptStat(0);

  string filename = "biasValues/pull_error_";
  filename += Form("PhoPt_%dto%d", (int) minpt, (int) maxpt);
  filename += "_" + cat;

  if (minpt >= 15. && maxpt <= 20.) filename += "_QCD30to50_Signal_Zg";
  else filename += "_QCD50to80_Signal_Zg";

  if (ele) filename += "_ele.txt";
  else filename += "_mu.txt";

  cout << filename << endl;

  ifstream infile;


  TH2F *hbias_sb = new TH2F("hbias_sb", "scan iso", 8, 3, 11, 11, 4, 15);
  TH2F *herr_sb = new TH2F("herr_sb", "error of template", 8, 3, 11, 11, 4, 15);

  //loop to fill histogram
  for (int ix = 1; ix <= 8; ix++) {
    for (int iy = ix; iy <= 11; iy+=2){

      infile.open(filename.c_str());

      if (!infile) {
	cout << "file could not be opened" << endl;
	return (0);
      }

      if (infile.is_open()) {
	string out;
	while (getline(infile, out)) {
	stringstream ss(out);
	
	string s1, s2, s3;
	int a1, a2;
	float bias, err;

	ss >> a1 >> s1 >> s2 >> s3 >> a2 >> bias >> err;
	//std::cout << a1 << "\t" << s1 << "\t " << s2 << "\t" << a2 << "\t" << s3 << "\t" << bias << endl;
	//std::cout << "bias of signal yield: "  << bias << endl;
	if ( a1 == hbias_sb->GetXaxis()->GetBinLowEdge(ix) && a2 == hbias_sb->GetYaxis()->GetBinLowEdge(iy)+1) {
	  hbias_sb->SetBinContent(ix, iy, bias);
	  herr_sb->SetBinContent(ix, iy, err);
	//std::cout << "lower cut: " << a1 << "\t bin x: " << hbias_sb->GetXaxis()->GetBinLowEdge(ix) << endl;
	}
	}
      }
      infile.close();
    }
  }


  //find minimum and maximum value
  float min = 99.;
  float max = -99.;

  infile.open(filename.c_str());
  if (!infile) {
    cout << "file could not be opened" << endl;
    return (0);
  }

  if (infile.is_open()) {
    string out;
    while (getline(infile, out)) {
      stringstream ss(out);
      
      string s1, s2, s3;
      int a1, a2;
      float bias;
      
      ss >> a1 >> s1 >> s2 >> s3 >> a2 >> bias;
      if (min > bias) min = bias;
      if (max < bias) max = bias;
    }
  }

  cout << "min bias is: " << min << endl;
  cout << "max bias is: " << max << endl;
  //if ( abs(max/min) > 20) max = max/2;
  //cout << "max: " << max << endl;

  //for drawing histogram

  TCanvas *c1 = new TCanvas("c1", "bias of signal", 650, 650);
  c1->cd();
  c1->SetRightMargin(0.13);
  hbias_sb->GetYaxis()->SetTitle("PFChbias Upper Cut [GeV]");
  hbias_sb->GetXaxis()->SetTitle("PFChbias Lower Cut [GeV]");
  //hbias_sb->GetZaxis()->SetRangeUser(min, max);
  hbias_sb->SetMarkerSize(1.5);
  hbias_sb->Draw("colz text");

  TCanvas *cerr = new TCanvas("cerr", "error of template", 650, 650);
  cerr->cd();
  cerr->SetRightMargin(0.13);
  herr_sb->GetYaxis()->SetTitle("PFCherr Upper Cut [GeV]");
  herr_sb->GetXaxis()->SetTitle("PFCherr Lower Cut [GeV]");
  //herr_sb->GetZaxis()->SetRangeUser(min, max);
  herr_sb->SetMarkerSize(1.5);
  herr_sb->Draw("colz text");


  //return 0;
}
