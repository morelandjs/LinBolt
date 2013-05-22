#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "scattering.h"

scattering::scattering(ParameterReader* _paraRdr)
{
   paraRdr = _paraRdr;
   alphas = paraRdr->getVal("alphas");
   gs = pow(4.0*M_PI*alphas,0.5);
   sm = paraRdr->getVal("sm"); 
   gq=6.0;
}

scattering::~scattering()
{
}

double scattering::Gamma_i(double E1, double T, int i)
{
  double G=0.0;
  for(int j=-3; j<=3; j++){
  for(int k=-3; k<=3; k++){
  for(int l=-3; l<=3; l++){
    G += Gamma_ij2kl(E1,T,i,j,k,l);
  }
  }
  }
  return G;
}

double scattering::Gamma_ij2kl(double E1, double T, int i, int j, int k, int l)
{
  double mreg = sm*gs*T;
  double G=0.0, dG=0.0;
  double E2min = pow(mreg,2.0)/(2.0*E1); 
  double E2max = 100; 
  int npart = 10000;
  double dE2 = (E2max-E2min)/npart;
  double E2 = E2min;
  for(int n=0; n<npart; n++){
    dG = dE2*F(E2,T)*Sigma_ij2kl(E1,E2,mreg,i,j,k,l);
    G += dG;
    E2 += dE2;
  }
  return G/(16.0*pow(M_PI,2.0)*pow(E1,2.0));
}

double scattering::Sigma_ij2kl(double E1, double E2, double m, int i, int j, int k, int l)
{
  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return 1.0/(32.0*M_PI)*(9.0/2.0)*pow(gs,4.0)*(3.0*A10(E1,E2,m)-A11(E1,E2,m)-A12(E1,E2,m)-A12(E1,E2,m));
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && l==-k){
    return 1.0/(16.0*M_PI)*(3.0/8.0)*pow(gs,4.0)*((4.0/9.0)*A14(E1,E2,m)-A13(E1,E2,m));
  }
  // {g q -> g q, g qbar -> g qbar}
  else if(i==0 && abs(j)>0 && k==0 && l==j){
    return 1.0/(16.0*M_PI)*pow(gs,4.0)*(A15(E1,E2,m)-(4.0/9.0)*A16(E1,E2,m));
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==k && j==l && i!=j && k!=l){
    return 1.0/(16.0*M_PI)*(4.0/9.0)*pow(gs,4.0)*A15(E1,E2,m);
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==j && j==k && k==l){
    return 1.0/(32.0*M_PI)*(4.0/9.0)*pow(gs,4.0)*(A15(E1,E2,m)+A15(E1,E2,m)-(2.0/3.0)*A17(E1,E2,m));
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i!=k && j!=l){
    return 1.0/(16.0*M_PI)*(4.0/9.0)*pow(gs,4.0)*A13(E1,E2,m);
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i==k && j==l){
    return 1.0/(16.0*M_PI)*(4.0/9.0)*pow(gs,4.0)*(A15(E1,E2,m)+A13(E1,E2,m)-(2.0/3.0)*A18(E1,E2,m));
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(j)>0 && k==0 && l==0 && i==-j){
    return 1.0/(32.0*M_PI)*(8.0/3.0)*pow(gs,4.0)*((4.0/9.0)*A14(E1,E2,m)-A13(E1,E2,m));
  }
  // >>no such event<<
  else{
    return 0;
  }
}

double scattering::F(double E2, double T)
{
  return gq*pow(exp(E2/T)+1,-1.0);
}

double scattering::A10(double E1, double E2, double m)
{
  return 4.0*E1*E2-2.0*pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))-2.0*pow(m,2.0);
}

double scattering::A11(double E1, double E2, double m)
{
  return (2.0/3.0)*E1*E2-(3.0/4.0)*pow(m,2)+pow(m,4)/(4.0*E1*E2)-pow(m,6)/(48.0*pow(E1,2)*pow(E2,2));
}

double scattering::A12(double E1, double E2, double m)
{
  return -8.0*pow(E1,2.0)*pow(E2,2.0)/pow(m,2.0)+4.0*E1*E2*log(4.0*E1*E2/pow(m,2.0)-1.0)+2.0*pow(m,2.0);
}

double scattering::A13(double E1, double E2, double m)
{
  return (8.0/3.0)*E1*E2-2.0*pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))-(1.0/2.0)*pow(m,2.0)-pow(m,4.0)/(2.0*E1*E2)+pow(m,6.0)/(24.0*pow(E1,2.0)*pow(E2,2.0));
}

double scattering::A14(double E1, double E2, double m)
{
  return 2.0*(4.0*E1*E2-pow(m,2.0))*log(4.0*E1*E2/pow(m,2.0)-1.0)-16.0*E1*E2+4.0*pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))+8.0*pow(m,2.0);
}

double scattering::A15(double E1, double E2, double m)
{
  return 16.0*pow(E1,2.0)*pow(E2,2.0)/pow(m,2.0)-8.0*E1*E2*log(4.0*E1*E2/pow(m,2.0)-1.0)+4.0*E1*E2-2.0*pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))-6.0*pow(m,2.0);
}

double scattering::A16(double E1, double E2, double m)
{
  return -(4.0*E1*E2-pow(m,2.0))*log(4.0*E1*E2/pow(m,2.0)-1.0)+2.0*E1*E2+pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))-pow(m,2.0);
}

double scattering::A17(double E1, double E2, double m)
{
  return 2.0*(4.0*E1*E2-pow(m,2.0))*log(4.0*E1*E2/pow(m,2.0)-1.0)-8.0*E1*E2+4.0*pow(m,2.0);
}

double scattering::A18(double E1, double E2, double m)
{
  return -(4.0*E1*E2-pow(m,2.0))*log(4.0*E1*E2/pow(m,2.0)-1.0)+10.0*E1*E2-3.0*pow(m,2.0)*log(2.0*E1*E2/pow(m,2.0))-5.0*pow(m,2.0);
}

void scattering::printGamma(char filename[])
{
  ofstream outfile;
  outfile.open(filename, std::ios_base::trunc);

  double Tmin = 0.1;
  double Tmax = 1.5; 
  int npart = 1000;
  double dT = (Tmax-Tmin)/npart;
  double T = Tmin;

  /*double mreg = sm*gs*0.5;
  double E2min = pow(mreg,2.0)/(2.0*5.0);
  double E2max = 5.0; 
  int npart = 1000;
  double dE2 = (E2max-E2min)/npart;
  double E2 = E2min;*/

  for(int n=0; n<npart; n++){
    outfile << setprecision(12) << setw(22) << T 
	    << setprecision(12) << setw(22) << Gamma_i(5.0,T,0) << endl;
    T += dT;
  }
  outfile.close();
}
