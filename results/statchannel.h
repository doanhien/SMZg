#ifndef statchannel_h
#define statchannel_h

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "TGraphErrors.h"

#include "lumierror.h"

using namespace std;

enum asym {
   XS_EB=0,
   XS_EE,
   XS_TOT,
};

class statchannel
{
   private:
      string m_label;
      int m_nbins;
      vector<double> m_bins;
      vector<double> m_eb_yields;
      vector<double> m_eb_staterr;
      vector<double> m_ee_yields;
      vector<double> m_ee_staterr;
      vector<double> m_eff;
      vector<double> m_efferr;
      //vector< vector<double> > m_eff;
      //vector< vector<double> > m_efferr;
      vector< vector<double> > m_eff_sf;
      vector<double> m_th;
      vector<double> m_therr;

   public:
      statchannel() {};
      ~statchannel() {};

      void read(int n, const char* filename);
      void print();
      
      // getters
      string get_label() const {return m_label;};
      int get_nbins() const {return m_nbins;};
      vector<double> get_bins() const {return m_bins;};
      vector<double> get_eb_yields() const {return m_eb_yields;};
      vector<double> get_eb_staterr() const {return m_eb_staterr;};
      vector<double> get_ee_yields() const {return m_ee_yields;};
      vector<double> get_ee_staterr() const {return m_ee_staterr;};
      vector<double> get_eff() const {return m_eff;};
      vector<double> get_efferr() const {return m_efferr;};
      vector< vector<double> > get_eff_sf() const {return m_eff_sf;};
      vector<double> get_th() const {return m_th;};
      vector<double> get_therr() const {return m_therr;};

      // static
      static TGraphErrors* graph(const statchannel chanp, asym mode=XS_EB, bool isth=true, bool dosyst=false);
      static double total(int n, vector<double> v);
      static double total2(int n, vector<double> v);
};

#endif // #ifdef statchannel_h
