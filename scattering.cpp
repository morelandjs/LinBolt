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

vector<int> scattering::sample2to2(double s, double temp, int i)
{
  vector<int> ij2kl(4,0);
  vector<double> partitions; partitions.push_back(0.0);
  double prob=0.0, total_prob=0.0, partition = 0.0; 
  double channel_rate=0.0, total_rate = Gamma_i(s,temp,i);
  double rand = drand48();

  cout << "total rate: " << total_rate << endl;
  cout << "rand: " << rand << endl;

  int index = 0;
  for(int j=-3; j<=3; j++){
    for(int k=-3; k<=3; k++){
      for(int l=-3; l<=3; l++){
	ij2kl[0] = i; ij2kl[1] = j; ij2kl[2] = k; ij2kl[3] = l;
	channel_rate = Gamma_ij2kl(s,temp,ij2kl);
	prob = channel_rate/total_rate;
	if(channel_rate > 0){
	  total_prob += prob;
	  cout << "total probability: " << total_prob << "  " << vector2string(ij2kl) << endl;
	  partition += prob;
	  partitions.push_back(partition);
	  if(partitions[index] < rand && rand < partitions[index+1]){
	    return ij2kl;
	  }
          index++;
	}
      }
    }
  }	
  cout << ">>error in decay channel sampling<<" << endl;
}

double scattering::Gamma_i(double E1, double temp, int i)
{
  int index=0;
  double G=0.0;
  vector<int> ij2kl(4,0);
  for(int j=-3; j<=3; j++){
    for(int k=-3; k<=3; k++){
      for(int l=-3; l<=3; l++){
	ij2kl[0]=i;
	ij2kl[1]=j;
	ij2kl[2]=k;
	ij2kl[3]=l;
	G += Gamma_ij2kl(E1,temp,ij2kl);
	index++;
      }
    }
  }
  return G;
}

double scattering::Gamma_ij2kl(double E1, double temp, vector<int> ij2kl)
{
  double mreg = sm*gs*temp;
  double G=0.0, dG=0.0;
  int npart = 10000;
  double E2min = pow(mreg,2.0)/(2.0*E1), E2max = 50.0; 
  double E2 = E2min, dE2 = (E2max-E2min)/npart;

  for(int n=0; n<npart; n++){
    dG = dE2*F(E2,temp)*max(Sigma_ij2kl(E1,E2,mreg,ij2kl),0.0);
    G += dG;
    E2 += dE2;
  }
  return G/(16.0*pow(M_PI,2.0)*pow(E1,2.0));
}

double scattering::Sigma_ij2kl(double E1, double E2, double m, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

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
    return -1.0;
  }
}

double scattering::sampleTheta(double s, double temp, vector<int> ij2kl)
{
  double mreg = sm*gs*temp;
  double tmin = -s+(mreg*mreg);
  double tmax = -mreg*mreg;
  double dIM2 = (IM2_ij2kl(s,tmax,ij2kl)-IM2_ij2kl(s,tmin,ij2kl));
  double dIM20 = drand48()*dIM2;
  double IM20 = (IM2_ij2kl(s,tmin,ij2kl)+dIM20);
  double t = invert_IM2_ij2kl(s,temp,ij2kl,IM20); 
  double theta = acos(2.0*t/s+1.0);
  return theta;
}

double scattering::samplePhi()
{
  double phi = drand48()*2.0*M_PI;
  return phi;
}

double scattering::invert_IM2_ij2kl(double s, double temp, vector<int> ij2kl, double IM2)
{
  double mreg = sm*gs*temp;
  double tmin = -s+(mreg*mreg);
  double tmax = -mreg*mreg;
  double dist = (IM2_ij2kl(s,tmax,ij2kl)-IM2_ij2kl(s,tmin,ij2kl));
  double precision = 0.1;
  double t = tmin+drand48()*(tmax-tmin);
  do{
    dist = IM2_ij2kl(s,t,ij2kl)-IM2;
    if(dist<0.0){
      tmin = t;
      t = t+drand48()*(tmax-t);
    }
    else{
      tmax = t;
      t = tmin+drand48()*(t-tmin);
    }
  }while(abs(dist) > precision); 
  return t;
}

double scattering::IM2_ij2kl(double s, double t, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

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
    return -1;
  }
}

double scattering::M2_ij2kl(double s, double t, double u, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

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
    return -1;
  }
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

double scattering::getMandelstamS(vector<double> &p1, vector<double> &p2)
{
  double v0, v1, v2, v3, s;
  v0 = pow(p1[0]+p2[0],2.0);
  v1 = pow(p1[1]+p2[1],2.0);
  v2 = pow(p1[2]+p2[2],2.0);
  v3 = pow(p1[3]+p2[3],2.0);

  s = v0-v1-v2-v3;
  return s;
}

double scattering::F(double E2, double temp)
{
  return gq*pow(exp(E2/temp)+1,-1.0);
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
