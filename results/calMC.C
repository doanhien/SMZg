#include "math.h"

using namespace std;

void calMC(const char* filename = "yield_mc_eb_ele.txt") {


  int m_nbins = 14;
  ifstream file;
  file.open(filename);
  double x;

  vector<double> m_bins;
  string m_label;
  vector<double> m_yields;
  vector<double> m_staterr;
  vector< vector<double> > m_eff;
  //vector< vector<double> > m_efferr;

  //getline(file,m_label);
  m_bins.resize(m_nbins+1);

  cout << "start program" << endl;
  for (int i=0; i<m_nbins; i++) file >> m_yields[i];
  //m_staterr.resize(m_nbins);
  cout << "get error" << endl;
  for (int i=0; i<m_nbins; i++) file >> m_staterr[i];

  file >> x;
  while (!file.eof())
    {
      vector<double> tmpvec;
      tmpvec.resize(m_nbins);
      for (int i=0; i<m_nbins; i++) {tmpvec[i] = x; file >> x;};
      m_eff.push_back(tmpvec);
    }



  cout << "yields: " << endl;
  for (int i=0; i<m_yields.size(); i++) cout << m_yields[i] << " ";
  cout << endl;
  cout << "statistical error on yields: " << endl;
  for (int i=0; i<m_staterr.size(); i++) cout << m_staterr[i] << " ";
  cout << endl;
  for (int j=0; j<m_eff.size(); j++)
    {
      cout << "eff. " << j << ": " << endl;
      for (int i=0; i<m_yields.size(); i++) cout << m_eff[j][i] << " ";
      cout << endl;
    }


  for (int i=0; i<m_yields.size(); i++) {
    for (int j=0; j<m_eff.size(); j++) {
      m_yields[i] /= m_eff[j][i];
      m_staterr[i] /= m_eff[j][i];
    }

    cout << "theory: " << m_yields[i] << endl;
    cout << "theory error: " << m_staterr[i] << endl;
  }

}


