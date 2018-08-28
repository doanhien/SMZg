#include "statchannel.h"

#include "math.h"

#include <iomanip>

using namespace std;

void statchannel::read(int n, const char* filename)
{
   m_nbins = n;
   ifstream file;
   file.open(filename);
   double x;

   getline(file,m_label);
   m_bins.resize(m_nbins+1);
   for (int i=0; i<m_nbins+1; i++) file >> m_bins[i];
   m_eb_yields.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_eb_yields[i];
   m_eb_staterr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_eb_staterr[i];
   m_ee_yields.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_ee_yields[i];
   m_ee_staterr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_ee_staterr[i];
   m_th.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_th[i];
   m_therr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_therr[i];
   m_eff.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_eff[i];
   m_efferr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_efferr[i];

   file >> x;
   while (!file.eof())
   {
      vector<double> tmpvec;
      tmpvec.resize(m_nbins);
      for (int i=0; i<m_nbins; i++) {tmpvec[i] = x; file >> x;};
      m_eff_sf.push_back(tmpvec);
      //for (int i=0; i<m_nbins; i++) {tmpvec[i] = x; file >> x;};
      //m_efferr.push_back(tmpvec);
   }
}

void statchannel::print()
{
   cout << m_label << endl;
   cout << "bin edges: " << endl;
   for (int i=0; i<m_nbins+1; i++) cout << m_bins[i] << " ";
   cout << endl;
   cout << "yields in EB: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_eb_yields[i] << " ";
   cout << endl;
   cout << "statistical error on yields in EB: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_eb_staterr[i] << " ";
   cout << endl;
   cout << "yields in EE: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_ee_yields[i] << " ";
   cout << endl;
   cout << "statistical error on yields in EE: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_ee_staterr[i] << " ";
   cout << endl;
   cout << "theory: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_th[i] << " ";
   cout << endl;
   cout << "theory error: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_therr[i] << " ";
   cout << endl;
   cout << "eff. " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_eff[i] << " ";
   cout << endl;
   cout << "error on eff. " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_efferr[i] << " ";
   cout << endl;
   for (int j=0; j<m_eff_sf.size(); j++)
   {
      cout << "eff with up and down SF. " << j << ": " << endl;
      for (int i=0; i<m_nbins; i++) cout << m_eff_sf[j][i] << " ";
      cout << endl;
   }
   cout << endl;
}

TGraphErrors* statchannel::graph(const statchannel chan, asym mode, bool isth, bool dosyst)
{

  //cout << "Start calculate xs and put to graph" << endl;
   int nbins = chan.get_nbins();
   vector<double> bins = chan.get_bins();
   vector<double> th = chan.get_th();
   vector<double> theb = chan.get_th();
   vector<double> thee = chan.get_th();

   vector<double> therr = chan.get_therr();
   vector<double> therr_sys = chan.get_therr();
   vector<double> theb_err = chan.get_therr();
   vector<double> thee_err = chan.get_therr();


   if (!isth)
   {
      theb = chan.get_eb_yields();
      theb_err = chan.get_eb_staterr();
      thee = chan.get_ee_yields();
      thee_err = chan.get_ee_staterr();
      vector<double > eff = chan.get_eff();
      vector<double > efferr = chan.get_efferr();
      vector< vector<double> > eff_sf = chan.get_eff_sf();

      vector<double > yield_sf;

      // now compute the errors and the yields taking into account the efficiencies
      bool correlm = mode!=XS_TOT;
      int neff = eff_sf.size();

      for (int i=0; i<nbins; i++)
      {
	th[i] = 0.,  therr[i] = 0, therr_sys[i] = 0;
	theb[i] /= eff[i];
	thee[i] /= eff[i];
	theb_err[i] /= eff[i];
	thee_err[i] /= eff[i];
	th[i] = theb[i] + thee[i] ;
	 
         if (!dosyst)
         {
	   therr[i] = pow(theb_err[i],2) + pow(thee_err[i],2);
	 }

         // now the systematics (be careful with correlations)
         if (dosyst)
         {
	   therr[i] = pow(theb_err[i],2) + pow(thee_err[i],2);
	   cout << "calculating systematics errors: " << endl;
            for (int ieff=0; ieff<neff-1; ieff+=2)
            {
	      float tmpmax = 0.;
	      float yield_sf_up;
	      float yield_sf_down;
	      yield_sf_up = th[i]*eff[i]/eff_sf[ieff][i];
	      yield_sf_down = th[i]*eff[i]/eff_sf[ieff+1][i];
	      if (fabs(th[i] - yield_sf_up) > fabs(th[i] - yield_sf_down)) tmpmax = fabs(th[i] - yield_sf_up);
	      else tmpmax = fabs(th[i] - yield_sf_down);
	      //cout << "efficiency up: " << eff_sf[ieff][i] << " efficiency down: " << eff_sf[ieff+1][i] << endl;
	      //cout << "yield_sf_up: " << yield_sf_up << " yield_sf_down: "<< yield_sf_down << " nominal yield: " << th[i] << endl;
	      therr_sys[i] += pow(tmpmax,2);
	      if (i < 4) cout << tmpmax*100/(th[i]) << " " ; 
	      else if ( i < 10) cout << setprecision(3)<< tmpmax*100/(th[i]) << " " ;
	      else if ( i < 11) cout << setprecision(3)<< tmpmax*100/(th[i]) << " " ;
	      else cout << setprecision(3) << tmpmax*100/(th[i]) << " " ;
            }
	    cout << "\n ---------------" << endl;

	 }

         // at last, take the square root of the sum of the squares.
	 therr_sys[i] += pow(th[i]*LUMIERR,2); //for lumi
         therr[i] = sqrt(therr[i]);
         therr_sys[i] = sqrt(therr_sys[i]);
	 //if (!isth && dosyst) cout << "individual systematics un: " << therr_sys[i] << endl;
      }
   }

   double *x, *y, *ex, *ey;
   x = new double[nbins];
   y = new double[nbins];
   ex = new double[nbins];
   ey = new double[nbins];

   for (int i=0; i<nbins; i++)
     //for (int i=0; i<5; i++) //for blinded
   {

     int ix = i;
     x[i] = (bins[ix+1]+bins[ix])/2.;
     ex[i] = (bins[ix+1]-bins[ix])/2.;

      switch (mode)
      {
         case XS_EB : // yields in EB only
            {
               y[i] = theb[i]/(bins[ix+1]-bins[ix]);
               ey[i] = theb_err[i]/(bins[ix+1]-bins[ix]);

               break;
            }

         case XS_EE : // yields in EE only
            {
               y[i] = thee[i]/(bins[ix+1]-bins[ix]);
               ey[i] = thee_err[i]/(bins[ix+1]-bins[ix]);
               break;
            }

         case XS_TOT : // xs of EB+EE
            {
               y[i] = th[i]/(bins[ix+1]-bins[ix]);
	       ey[i] = therr[i];
	       if ( !isth && !dosyst) cout << "statistical error:" << ey[i] << endl;

               // for the systematic error we need to be a bit careful with correlations
	       if (!isth && dosyst)
               {
		 ey[i] = sqrt(pow(therr[i],2) + pow(therr_sys[i],2));
		 //cout << "systematics error:" << therr_sys[i]/((bins[ix+1]-bins[ix])*35.9) << endl;
		 //cout << "systematics error:" << therr_sys[i]*100/th[i] << endl;
		 cout << "stat error: " << therr[i] << "\t systematics error:" << therr_sys[i] << endl;
               }
               ey[i] = ey[i]/(bins[ix+1]-bins[ix]);
               break;
            }

         default :
            {
               cout << "statchannel::graph_th: Error, unknown asymmetry type!"  << endl;
               return NULL;
            }
      }
      if (!isth) {
	y[i] /= LUMI;
	ey[i] /= LUMI;
      }
   }


   //cout << "put into the graph " << endl;   
   TGraphErrors *gr = new TGraphErrors(nbins,x,y,ex,ey);
   if (mode == XS_TOT) gr->Print();
   cout << "end !!!" << endl;
   return gr;
}

double statchannel::total(int n, vector<double> v)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += v[i];
   return tot;
}

double statchannel::total2(int n, vector<double> v)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += v[i]*v[i];
   return tot;
}
