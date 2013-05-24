#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
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

double scattering::sampleDiffXS(double s, double temp, int i, int j, int k, int l)
{
  double mreg = sm*gs*temp;
  double tmin = -s+(mreg*mreg);
  double tmax = -mreg*mreg;
  double IM2min = 0.0;
  double IM2max = (M2_ij2kl_integrated(s,tmax,i,j,k,l)-M2_ij2kl_integrated(s,tmin,i,j,k,l));
  double IM2rand = drand48()*IM2max;
  // find t for given sigma value
  double IM2targ = (M2_ij2kl_integrated(s,tmin,i,j,k,l)+IM2rand);
  double t = invertDiffXS(s,temp,i,j,k,l,IM2targ);
  return acos(2.0*t/s+1.0);
}

double scattering::invertDiffXS(double s, double temp, int i, int j, int k, int l, double IM2targ)
{
  double mreg = sm*gs*temp;
  double tmin = -s+(mreg*mreg);
  double tmax = -mreg*mreg;
  double dist = (M2_ij2kl_integrated(s,tmax,i,j,k,l)-M2_ij2kl_integrated(s,tmin,i,j,k,l));
  double precision = 0.00001;
  double trand = tmin+drand48()*(tmax-tmin);
  double tlb = tmin, tub = tmax;
  do{
    dist = M2_ij2kl_integrated(s,trand,i,j,k,l)-IM2targ;
    if(dist<0.0){
      tlb = trand;
      trand = trand+drand48()*(tub-trand);
    }
    else{
      tub = trand;
      trand = tlb+drand48()*(trand-tlb);
    }
  }while(abs(dist) > precision); 
  return trand;
}

double scattering::Gamma_i(double E1, double temp, int i)
{
  double G=0.0;
  for(int j=-3; j<=3; j++){
  for(int k=-3; k<=3; k++){
  for(int l=-3; l<=3; l++){
    G += Gamma_ij2kl(E1,temp,i,j,k,l);
  }
  }
  }
  return G;
}

double scattering::Gamma_ij2kl(double E1, double temp, int i, int j, int k, int l)
{
  double mreg = sm*gs*temp;
  double G=0.0, dG=0.0;
  double E2min = pow(mreg,2.0)/(2.0*E1); 
  double E2max = 100; 
  int npart = 10000;
  double dE2 = (E2max-E2min)/npart;
  double E2 = E2min;
  for(int n=0; n<npart; n++){
    dG = dE2*F(E2,temp)*Sigma_ij2kl(E1,E2,mreg,i,j,k,l);
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

double scattering::M2_ij2kl_integrated(double s, double t, int i, int j, int k, int l)
{
  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return (9.0*pow(gs,4.0)*(-(pow(s,4.0)/t) + 3.0*pow(s,2.0)*t + (s*pow(t,2.0))/2.0 + pow(t,3.0)/3.0 - pow(s,4.0)/(s+t) + pow(s,3.0)*log(abs(t)) - pow(s,3.0)*log(abs(s+t))))/(2.0*pow(s,2.0));
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && l==-k){
    return -((pow(gs,4.0)*(17.0*pow(s,2.0)*t + 9.0*s*pow(t,2.0) + 6.0*pow(t,3.0) + 4.0*pow(s,3.0)*log(abs(t)) - 4.0*pow(s,3.0)*log(abs(s+t))))/(24.0*pow(s,2.0)));
  }
  // {g q -> g q, g qbar -> g qbar}
  else if(i==0 && abs(j)>0 && k==0 && l==j){
    return (pow(gs,4.0)*(-((18.0*pow(s,3.0))/t) + 13.0*s*t + 2.0*pow(t,2.0) + 18.0*pow(s,2.0)*log(abs(t)) + 4.0*pow(s,2.0)*log(abs(s+t))))/(9.0*s);
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==k && j==l && i!=j && k!=l){
    return 4.0/9.0*pow(gs,4.0)*(-((2.0*pow(s,2.0))/t) + t + 2*s*log(abs(t)));
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==j && j==k && k==l){
    return 8.0/27.0*pow(gs,4.0)*(-((3.0*pow(s,2.0))/t) + 3.0*t - (3.0*pow(s,2.0))/(s+t) + 4.0*s*log(abs(t)) - 4.0*s*log(abs(s+t)));
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i!=k && j!=l){
    return (4.0*pow(gs,4.0)*(pow(s,2.0)*t + s*pow(t,2.0) + (2.0*pow(t,3.0))/3.0))/(9.0*pow(s,2.0));
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i==k && j==l){
    return (8.0*pow(gs,4.0)*(-((3.0*pow(s,4.0))/t) + pow(s,2.0)*t + s*pow(t,2.0) + pow(t,3.0) + 2.0*pow(s,3.0)*log(abs(t))))/(27.0*pow(s,2.0));
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(j)>0 && k==0 && l==0 && i==-j){
    return -((8.0*pow(gs,4.0)*(17.0*pow(s,2.0)*t + 9.0*s*pow(t,2.0) + 6.0*pow(t,3.0) + 4.0*pow(s,3.0)*log(abs(t)) - 4.0*pow(s,3.0)*log(abs(s+t))))/(27.0*pow(s,2.0)));
  }
  // >>no such event<<
  else{
    return 0;
  }
}

double scattering::M2_ij2kl(double s, double t, double u, int i, int j, int k, int l)
{
  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return (9.0/2.0)*pow(gs,4.0)*(3.0-u*t/(s*s)-(u*s)/(t*t)-(s*t)/(u*u));
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && l==-k){
    return (3.0/8.0)*pow(gs,4.0)*((4.0/9.0)*(t*t+u*u)/(t*u)-(t*t+u*u)/(s*s));
  }
  // {g q -> g q, g qbar -> g qbar}
  else if(i==0 && abs(j)>0 && k==0 && l==j){
    return pow(gs,4.0)*((s*s+u*u)/(t*t)-(4.0/9.0)*(s*s+u*u)/(s*u));
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==k && j==l && i!=j && k!=l){
    return (4.0/9.0)*pow(gs,4.0)*(s*s+u*u)/(t*t);
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==j && j==k && k==l){
    return (4.0/9.0)*pow(gs,4.0)*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u)-(2.0/3.0)*(s*s)/(t*u));
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i!=k && j!=l){
    return (4.0/9.0)*pow(gs,4.0)*(t*t+u*u)/(s*s);
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(j)>0 && abs(k)>0 && abs(l)>0 && i==-j && k==-l && i==k && j==l){
    return (4.0/9.0)*pow(gs,4.0)*((s*s+u*u)/(t*t)+(t*t+u*u)/(s*s)-(2.0/3.0)*(u*u)/(s*t));
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(j)>0 && k==0 && l==0 && i==-j){
    return (8.0/3.0)*pow(gs,4.0)*((4.0/9.0)*(t*t+u*u)/(t*u)-(t*t+u*u)/(s*s));
  }
  // >>no such event<<
  else{
    return 0;
  }
}

double scattering::F(double E2, double temp)
{
  return gq*pow(exp(E2/temp)+1,-1.0);
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
  double temp = Tmin;

  /*double mreg = sm*gs*0.5;
  double E2min = pow(mreg,2.0)/(2.0*5.0);
  double E2max = 5.0; 
  int npart = 1000;
  double dE2 = (E2max-E2min)/npart;
  double E2 = E2min;*/

  for(int n=0; n<npart; n++){
    outfile << setprecision(12) << setw(22) << temp 
	    << setprecision(12) << setw(22) << Gamma_i(5.0,temp,0) << endl;
    temp += dT;
  }
  outfile.close();
}
