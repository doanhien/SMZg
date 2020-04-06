#include <iostream>
#include <fstream>

using namespace std;

void transferColRow(TString cat = "pt", string filename = "results/signal_yield_Ptg_EE_SB_6to14_ele_newZG.txt") {

  //string filename = "results/signal_yield_Ptg_EE_SB_6to14_ele_newZG.txt";
  //string filename = "results/signal_yield_Mllg_EB_newZG_newSB_ele_TUnfold_bin_isPVGood_biasUp.txt";

  cout << "iput file: " << filename << endl;
  ifstream fin;
  fin.open(filename.c_str());

  if (!fin) {
    cout << "file could not be opened" << endl;
    //return (0);
  }

  string out;

  cout << setprecision(2) << std::fixed << endl;
  cout << "yield: " << endl;
  while (getline(fin, out)) {
    //if (out.find("pt") != string::npos) {
    if (out.find(cat) != string::npos) {
      stringstream ss(out);
      string s1, s2, s3, s4, s5, s6, s7, s8, s9;
      int a1, a2;
      float yield, err, corr, err_corr, fk, err_fk;

      //ss >> a1 >> s1 >> s2  >> s3 >> a2 >> yield >> s4 >> err;
      ss >> a1 >> s1 >> s2  >> s3 >> a2 >> yield >> s4 >> err >> s5 >> s6 >> corr >> s7 >> err_corr >>s8 >> fk;
      //ss >> a1 >> s1 >> s2  >> s3 >> a2 >> yield >> s4 >> err >> s5 >> s6 >> corr >> s7 >> err_corr;
      cout << corr << " ";
    }
  }

  cout << endl;
  fin.close();

  fin.open(filename.c_str());
  cout << "error of yield: " << endl;
  while (getline(fin, out)) {
    if (out.find(cat) != string::npos) {
      stringstream ss(out);
      string s1, s2, s3, s4, s5, s6, s7,s8, s9;
      int a1, a2;
      float yield, err, corr, err_corr, fk, err_fk;
      ss >> a1 >> s1 >> s2  >> s3 >> a2 >> yield >> s4 >> err >> s5 >> s6 >> corr >> s7 >> err_corr >> s8 >> fk >> s9 >>err_fk;
      //ss >> a1 >> s1 >> s2  >> s3 >> a2 >> yield >> s4 >> err >> s5 >> s6 >> corr >> s7 >> err_corr;
      cout << err_corr << " ";
    }
  }

  cout << endl;

  cout << "done!!!" << endl;


}
