#include <iostream>     // std::cout
#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int biasCal (bool ele = true, string cat = "EB", float minpt = 15., float maxpt = 20.) {

  string filename = "biasValues/pull_error_";
  filename += Form("PhoPt_%dto%d", (int) minpt, (int) maxpt);
  filename += "_" + cat;

  if (minpt >= 15. && maxpt < =20.) filename += "_QCD30to50_Signal_Zg";
  else filename += "_QCD50to80_Signal_Zg";

  if (ele) filename += "_ele.txt";
  else filename += "_mu.txt";

  cout << filename << endl;

  ifstream infile;
  infile.open(filename.c_str());

  if (!infile) {
    cout << "file could not be opened" << endl;
    return (0);
    //continue;
  }

  string check_st1 = "nominal:" ;
  string check_st2 = "sig:" ;
  float sig_nom;

  vector<float> bias_lower_iso3;
  vector<float> bias_lower_iso4;
  vector<float> bias_lower_iso5;
  vector<float> bias_lower_iso6;
  vector<float> bias_lower_iso7;
  vector<float> bias_lower_iso8;
  vector<float> bias_lower_iso9;
  vector<float> bias_lower_iso10;


  TH2F *hiso_sb = new TH2F("hiso_sb", "scan iso", 12, 3, 15, 10, 5, 15);

  //loop to fill histogram
  for (int ix = 1; ix <= 12; ix++) {
    for (int iy = ix; iy <= 10; iy++){

      if (infile.is_open()) {
	string out;
	while (getline(infile, out)) {
	  //if (out.find(check_st1) != string::npos) {
	//cout << "find nomial yield it!" << endl;
	stringstream ss(out);
	
	string s1, s2, s3;
	int a1, a2;
	ss >> a1 >> s1 >> s2 >> a2 >> s3 >> bias;
	//std::cout << "string: " << a1 << "\t" << a2 << "\t " << a3 << "\t" << a4 << "\t" << a5 << "\t" << a6 << "\t" << sig_nom << endl;
	std::cout << "bias of signal yield: "  << bias << endl;
	if ( a1 == hiso_sb->GetXaxis()->GetLowerBin(ix)) hiso_sb->Fill(ix, iy, bias);
      }

      if (out.find(check_st2) != string::npos) {
	//cout << "find yield from side-band bkg template!" << endl;
	stringstream ss(out);
	
	float d;
	string b1, b2, b3, b4, b5, b6, b7;
	ss >> b1 >> b2 >> b3 >> b4 >> b5 >> b6 >> b7 >> d;
	//std::cout << "nomial signal yield: "  << sig_nom << ",  signal yield from sb bkg  " 
	//	  << b1 << b2 << b3 << b4 << b5 << "\t" << d << endl;

	float bias = (d - sig_nom) / sig_nom *100;
	cout << "bias of " << b1 << b2 << b3 << b4 << b5 << ":\t"<< fixed << setprecision(1) << bias << " %" << endl;


	if (b1.find("3") != string::npos) bias_lower_iso3.push_back(bias);	
	else if (b1.find("4") != string::npos) bias_lower_iso4.push_back(bias);
	else if (b1.find("5") != string::npos) bias_lower_iso5.push_back(bias);
	else if (b1.find("6") != string::npos) bias_lower_iso6.push_back(bias);
	else if (b1.find("7") != string::npos) bias_lower_iso7.push_back(bias);
	else if (b1.find("8") != string::npos) bias_lower_iso8.push_back(bias);
	else if (b1.find("9") != string::npos) bias_lower_iso9.push_back(bias);
	else if (b1.find("10") != string::npos) bias_lower_iso10.push_back(bias);


      }
    }
  }


  //for drawing histogram
  const int nbinX = 8;
  const int nbinY = 11;
  float bias[nbinX][nbinY];

  float iso_lower[nbinX+1] = {3, 3.99, 4.99, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99};
  float iso_upper[nbinY+1] = {5, 5.99, 6.99, 7.99, 8.99, 9.99, 10.99, 11.99, 12.99, 13.99, 14.99, 15.99};

  TH2F *hbias = new TH2F("hbias", "signal yield bias", nbinX, iso_lower, nbinY, iso_upper);

  int n = 0;
  for (int j = 1; j <= nbinY; j+=2) {
    hbias->SetBinContent(1, j, bias_lower_iso3[j-1-n]);
    hbias->SetBinContent(2, j, bias_lower_iso4[j-1-n]);
    hbias->SetBinContent(3, j, bias_lower_iso5[j-1-n]);
    hbias->SetBinContent(4, j, bias_lower_iso6[j-1-n]);
    hbias->SetBinContent(5, j, bias_lower_iso7[j-1-n]);
    hbias->SetBinContent(6, j, bias_lower_iso8[j-1-n]);
    hbias->SetBinContent(7, j, bias_lower_iso9[j-1-n]);
    hbias->SetBinContent(8, j, bias_lower_iso10[j-1-n]);

    n++;
  }


  cout << "=================================================" << endl;
  cout << "         lower iso (GEV)         " << endl;
  cout << "|\t |";
  for (int i = 3; i <= 10; i++) {
    cout << " " << fixed << setprecision(1) << (float) i << " |\t";
  }
  cout << endl;
  //cout << "|\t | 3.0 |\t 4.0 |\t 5.0 |\t 6.0 |\t 7.0 |\t 8.0 |\t 9.0 |\t 10.0 |" << endl;
  for (int i = 5; i <= 15; i++) {
    cout << "| " << fixed << setprecision(1) << (float) i << "\t | " << endl;
  }

  TCanvas *c1 = new TCanvas();
  c1->cd();
  hbias->Draw("text");


  return 0;
}
