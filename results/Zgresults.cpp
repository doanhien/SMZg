#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TGraphErrors.h"

// #include "stattest_gaus.h"
#include "statchannel.h"

using namespace std;

int main(int argc, const char** argv)
{
   if (argc!=4 && argc!=5)
   {
      cout << "Usage:" << endl;
      cout << argv[0] << " n_etabins1 channel1_plus.txt channel1_minus.txt n_etabins2 channel2_plus.txt chennael2_minus.txt" << endl;
      return -1;
   }

   bool singlechan = argc==4;

   statchannel chan1, chan2;
   int nbins1, nbins2;
   int nsyst1, nsyst2;

   TFile *f = new TFile("graph.root","RECREATE");
   f->cd();

   // channel 1
   cout << "running 1 channel" << endl;
   nbins1 = atoi(argv[1]);
   chan1.read(nbins1,argv[2]);
   chan1.print();

   // create the graphs
   TGraphErrors *gyields_th_1 = statchannel::graph(chan1, XS_TOT, true); gyields_th_1->SetName("gyields_th_1");
   TGraphErrors *gyields_exp_1 = statchannel::graph(chan1, XS_TOT, false, true);  gyields_exp_1->SetName("gyields_exp_1"); 
   TGraphErrors *gyields_exp_statonly_1 = statchannel::graph(chan1, XS_TOT, false, false);  gyields_exp_statonly_1->SetName("gyields_exp_statonly_1"); 

   //TGraphErrors *gyieldsee_th_1 = statchannel::graph(chan1eb, chan1ee, XS_EE, true);  gyieldsee_th_1->SetName("gyieldsee_th_1");
   //TGraphErrors *gyieldsee_exp_1 = statchannel::graph(chan1eb, chan1ee, XS_EE, false, true);  gyieldsee_exp_1->SetName("gyieldsee_exp_1"); 
   //TGraphErrors *gyieldsee_exp_statonly_1 = statchannel::graph(chan1eb, chan1ee, XS_EE, false, false);  gyieldsee_exp_statonly_1->SetName("gyieldsee_exp_statonly_1"); 


   if (argc == 4 && string(argv[3]) == "muon") {
     gyields_th_1->SetName("gyields_th_2");
     gyields_exp_1->SetName("gyields_exp_2");
     gyields_exp_statonly_1->SetName("gyields_exp_statonly_2");
     }

   gyields_th_1->Write();
   gyields_exp_1->Write();
   gyields_exp_statonly_1->Write();

   //gyieldseb_exp_1->Print();


   // channel 2
   if (!singlechan)
   {
      nbins2 = atoi(argv[3]);
      chan2.read(nbins2,argv[4]);
      // chan2.print();

      // create the graphs
      TGraphErrors *gyields_th_2 = statchannel::graph(chan2, XS_TOT, true); gyields_th_2->SetName("gyields_th_2"); gyields_th_2->Write();
      TGraphErrors *gyields_exp_2 = statchannel::graph(chan2, XS_TOT, false, true);  gyields_exp_2->SetName("gyields_exp_2"); gyields_exp_2->Write();
      TGraphErrors *gyields_exp_statonly_2 = statchannel::graph(chan2, XS_TOT, false, false);  gyields_exp_statonly_2->SetName("gyields_exp_statonly_2"); gyields_exp_statonly_2->Write();
   }

   // // yields per channel
   // // make the stattest objects
   // stattest_gaus statg1p(nbins1,nsyst1), statg1m(nbins1,nsyst1);
   // stattest_gaus statg2p(nbins2,nsyst2), statg2m(nbins2,nsyst2);
   // statg1p.set_data(chan1eb.get_yields());
   // statg1p.set_data_stat(chan1eb.get_staterr());
   // statg1p.set_pred(chan1eb.get_th());
   // statg1m.set_data(chan1ee.get_yields());
   // statg1m.set_data_stat(chan1ee.get_staterr());
   // statg1m.set_pred(chan1ee.get_th());

   f->Write();
   f->Close();
}

// NB
// objets stat à faire 
// (corréler les syst selon leur signe), ne pas oublier la lumi
// -> 1 par asymétrie
// -> 1 par canal (1p, 1m, 2p, 2m)
// -> 1p+1m, 2p+2m
// -> 1p+1m+2p+2m
//
// automatiser un peu tout ca
// -> une classe par canal
// -> une fonction qui prend l'objet stat, calcule le chi2 et fait la distrib de PE
//
// Creation des graphes
// - fct de statchannel
// - une fct par type de graphe (yields, C, A1p, A1m, A2, A3, A4)
// - retourne le graphe
