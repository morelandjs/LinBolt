#include <fstream>
#include <iomanip> 
#include <stdio.h>
#include <cmath>
#include <vector>
#include <iostream>
#include "checks.h"
using std::vector;
using namespace std;

checks::checks(ParameterReader* _paraRdr)
{
  paraRdr = _paraRdr;
}

void checks::printGamma(int species)
{
   // gluons
  vector<int> gggg(4,0), ggqqbar(4,0), gqgq(4,0);
  gggg[0] = 0; gggg[1] = 0; gggg[2] = 0; gggg[3] = 0;
  ggqqbar[0] = 0; ggqqbar[1] = 0; ggqqbar[2] = 1; ggqqbar[3] = -1;
  gqgq[0] = 0; gqgq[1] = 1; gqgq[2] = 0; gqgq[3] = 1;

  // quarks
  vector<int> qiqjqiqj(4,0), qiqiqiqi(4,0), qiqibarqjqjbar(4,0), qiqibarqiqibar(4,0), qqbargg(4,0); 
  qiqjqiqj[0]=1; qiqjqiqj[1]=2; qiqjqiqj[2]=1; qiqjqiqj[3]=2;
  qiqiqiqi[0]=1; qiqiqiqi[1]=1; qiqiqiqi[2]=1; qiqiqiqi[3]=1;
  qiqibarqjqjbar[0]=1; qiqibarqjqjbar[1]=-1; qiqibarqjqjbar[2]=2; qiqibarqjqjbar[3]=-2;
  qiqibarqiqibar[0]=1; qiqibarqiqibar[1]=-1; qiqibarqiqibar[2]=1; qiqibarqiqibar[3]=-1;
  qqbargg[0]=1; qqbargg[1]=-1; qqbargg[2]=0; qqbargg[3]=0;

  // charm
  vector<int> cqcq(4,0), cgcg(4,0);
  cgcg[0]=4; cgcg[1]=0; cgcg[2]=4; cgcg[3]=0;
  cqcq[0]=4; cqcq[1]=1; cqcq[2]=4; cqcq[3]=1;

  int nE1 = 1, nT = 100;
  double E1_min = 5.0, E1_max = 5.0, dE1 = (E1_max-E1_min)/double(nE1);
  double T_min = 0.1, T_max = 1.1, dT = (T_max-T_min)/double(nT);
  
  double incl_rate, gggg_rate, ggqqbar_rate, gqgq_rate, qgqg_rate, qiqjqiqj_rate, 
    qiqiqiqi_rate, qiqibarqjqjbar_rate, qiqibarqiqibar_rate, qqbargg_rate,
    cgcg_rate, cqcq_rate; 
  
  medium* oscar = new medium(paraRdr, "/run/media/morelandjs/b09c7850-22ca-4fc5-b7a3-96981704d524/jsm55/OSCAR2008H.dat");
  scattering* dynamics = new scattering(paraRdr);

  // open output file
  ofstream Gamma_vs_T, Gamma_vs_E1;
  Gamma_vs_T.open("/home/morelandjs/Research/LinBolt/results/Gamma_vs_T.dat", std::ios_base::trunc);
  Gamma_vs_E1.open("/home/morelandjs/Research/LinBolt/results/Gamma_vs_E1.dat", std::ios_base::trunc);

  // vs T
  for(int iT=0; iT<nT; iT++)
    {
      double T = T_min + iT*dT;
      double E1 = E1_min;

      if(species == 0){
	gggg_rate = dynamics->Gamma(E1,T,gggg);
	gqgq_rate = 6.0*dynamics->Gamma(E1,T,gqgq);
	ggqqbar_rate = dynamics->Gamma(E1,T,ggqqbar);
	incl_rate = dynamics->incl_Gamma(E1,T,0);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << gggg_rate << setprecision(12) << setw(22) << gqgq_rate << setprecision(12) << setw(22) << ggqqbar_rate
		   << setprecision(12) << setw(22) << incl_rate << endl;
      }

      if(species == 1){
	qgqg_rate = dynamics->Gamma(E1,T,gqgq);
	qiqjqiqj_rate = 4.0*dynamics->Gamma(E1,T,qiqjqiqj);
	qiqiqiqi_rate = dynamics->Gamma(E1,T,qiqiqiqi);
	qiqibarqiqibar_rate = dynamics->Gamma(E1,T,qiqibarqiqibar);
	qiqibarqjqjbar_rate = dynamics->Gamma(E1,T,qiqibarqjqjbar);
	qqbargg_rate = dynamics->Gamma(E1,T,qqbargg);
	incl_rate = dynamics->incl_Gamma(E1,T,1);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << qgqg_rate << setprecision(12) << setw(22) << qiqjqiqj_rate << setprecision(12) << setw(22) << qiqiqiqi_rate
		   << setprecision(12) << setw(22) << qiqibarqiqibar_rate << setprecision(12) << setw(22) << qiqibarqjqjbar_rate << setprecision(12) << setw(22) << qqbargg_rate 
		   << setprecision(12) << setw(22) << incl_rate << endl;
      }

      if(species==4){
	cqcq_rate = 6.0*dynamics->Gamma(E1,T,cqcq);
	cgcg_rate = dynamics->Gamma(E1,T,cgcg);
	incl_rate = dynamics->incl_Gamma(E1,T,4);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << cqcq_rate << setprecision(12) << setw(22) << cgcg_rate << setprecision(12) << setw(22) << incl_rate << endl;
      }
    }

  // vs E1
  for(int iE1=0; iE1<nE1; iE1++)
    {
      double T = T_min;
      double E1 = E1_min + iE1*dE1;
      
      if(species == 0){
	gggg_rate = dynamics->Gamma(E1,T,gggg);
	gqgq_rate = 6.0*dynamics->Gamma(E1,T,gqgq);
	ggqqbar_rate = dynamics->Gamma(E1,T,ggqqbar);
	incl_rate = dynamics->incl_Gamma(E1,T,0);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << gggg_rate << setprecision(12) << setw(22) << gqgq_rate << setprecision(12) << setw(22) << ggqqbar_rate
		   << setprecision(12) << setw(22) << incl_rate << endl;
      }

      if(species == 1){
	qgqg_rate = dynamics->Gamma(E1,T,gqgq);
	qiqjqiqj_rate = 4.0*dynamics->Gamma(E1,T,qiqjqiqj);
	qiqiqiqi_rate = dynamics->Gamma(E1,T,qiqiqiqi);
	qiqibarqiqibar_rate = dynamics->Gamma(E1,T,qiqibarqiqibar);
	qiqibarqjqjbar_rate = dynamics->Gamma(E1,T,qiqibarqjqjbar);
	qqbargg_rate = dynamics->Gamma(E1,T,qqbargg);
	incl_rate = dynamics->incl_Gamma(E1,T,1);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << qgqg_rate << setprecision(12) << setw(22) << qiqjqiqj_rate << setprecision(12) << setw(22) << qiqiqiqi_rate
		   << setprecision(12) << setw(22) << qiqibarqiqibar_rate << setprecision(12) << setw(22) << qiqibarqjqjbar_rate << setprecision(12) << setw(22) << qqbargg_rate 
		   << setprecision(12) << setw(22) << incl_rate << endl;
      }

      if(species==4){
	cqcq_rate = 6.0*dynamics->Gamma(E1,T,cqcq);
	cgcg_rate = dynamics->Gamma(E1,T,cgcg);
	incl_rate = dynamics->incl_Gamma(E1,T,4);
	Gamma_vs_T << setprecision(12) << setw(22) << T << setprecision(12) << setw(22) << cqcq_rate << setprecision(12) << setw(22) << cgcg_rate << setprecision(12) << setw(22) << incl_rate << endl;
      }
    }
 
  Gamma_vs_T.close();
  Gamma_vs_E1.close();
}
