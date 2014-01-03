#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdlib.h> 
#include "scattering.h"
#include "checks.h"
#include <vector>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_expint.h> 
using std::vector;
using namespace std;


// Particle constructor and deconstructor
//################################################

scattering::scattering(ParameterReader* _paraRdr)
{
  paraRdr = _paraRdr;
  gluon = paraRdr->getVal("gluon");
  quark = paraRdr->getVal("quark");
  hquark = paraRdr->getVal("hquark");
  verbose = paraRdr->getVal("verbose");
  mc = paraRdr->getVal("charm_mass"); 
  gg = paraRdr->getVal("gluon_dof"); 
  gq = paraRdr->getVal("quark_dof"); 
  alphas = paraRdr->getVal("alphas");
  gs = pow(4.0*M_PI*alphas,0.5);

  nE1 = paraRdr->getVal("nE1");
  E1_min = paraRdr->getVal("E1_min");
  E1_max = paraRdr->getVal("E1_max");
  if(E1_min==E1_max) dE1 = 1.0; 
  else dE1 = (E1_max-E1_min)/double(nE1);

  nE2 = paraRdr->getVal("nE2");
  E2_min = paraRdr->getVal("E2_min");
  E2_max = paraRdr->getVal("E2_max");
  if(E2_min==E2_max) dE2 = 1.0; 
  else dE2 = (E2_max-E2_min)/double(nE2);

  nT = paraRdr->getVal("nT");
  T_min = paraRdr->getVal("T_min");
  T_max = paraRdr->getVal("T_max");
  if(T_min==T_max) dT = 1.0;
  else dT = (T_max-T_min)/double(nT);

  ns_table = paraRdr->getVal("ns");

  cout << endl << endl << "generating lookup tables" << endl << endl;

  // gluons
  vector<int> gggg(4,0), ggqqbar(4,0), gqgq(4,0);
  gggg[0] = 0; gggg[1] = 0; gggg[2] = 0; gggg[3] = 0;
  ggqqbar[0] = 0; ggqqbar[1] = 0; ggqqbar[2] = 1; ggqqbar[3] = -1;
  gqgq[0] = 0; gqgq[1] = 1; gqgq[2] = 0; gqgq[3] = 1;

  // quarks
  vector<int> qiqjqiqj(4,0), qiqiqiqi(4,0), qiqibarqjqjbar(4,0), qiqibarqiqibar(4,0), qqbargg(4,0), qgqg(4,0); 
  qiqjqiqj[0]=1; qiqjqiqj[1]=2; qiqjqiqj[2]=1; qiqjqiqj[3]=2;
  qiqiqiqi[0]=1; qiqiqiqi[1]=1; qiqiqiqi[2]=1; qiqiqiqi[3]=1;
  qiqibarqjqjbar[0]=1; qiqibarqjqjbar[1]=-1; qiqibarqjqjbar[2]=2; qiqibarqjqjbar[3]=-2;
  qiqibarqiqibar[0]=1; qiqibarqiqibar[1]=-1; qiqibarqiqibar[2]=1; qiqibarqiqibar[3]=-1;
  qqbargg[0]=1; qqbargg[1]=-1; qqbargg[2]=0; qqbargg[3]=0;
  qgqg[0] = 1; qgqg[1] = 0; qgqg[2] = 1; qgqg[3] = 0;

  // charm
  vector<int> cqcq(4,0), cgcg(4,0);
  cgcg[0]=4; cgcg[1]=0; cgcg[2]=4; cgcg[3]=0;
  cqcq[0]=4; cqcq[1]=1; cqcq[2]=4; cqcq[3]=1;
 
  // gluon tables
  if(gluon){
    Gamma_gggg = populate_Gamma(gggg); Omega_gggg = populate_OmegaI(gggg);
    Gamma_ggqqbar = populate_Gamma(ggqqbar); Omega_ggqqbar = populate_OmegaI(ggqqbar);
    Gamma_gqgq = populate_Gamma(gqgq); Omega_gqgq = populate_OmegaI(gqgq);
  }
  // light quark tables
  if(quark){
    Gamma_qiqjqiqj = populate_Gamma(qiqjqiqj); Omega_qiqjqiqj = populate_OmegaI(qiqjqiqj);
    Gamma_qiqiqiqi = populate_Gamma(qiqiqiqi); Omega_qiqiqiqi = populate_OmegaI(qiqiqiqi);
    Gamma_qiqibarqjqjbar = populate_Gamma(qiqibarqjqjbar); Omega_qiqibarqjqjbar = populate_OmegaI(qiqibarqjqjbar);
    Gamma_qiqibarqiqibar = populate_Gamma(qiqibarqiqibar); Omega_qiqibarqiqibar = populate_OmegaI(qiqibarqiqibar);
    Gamma_qqbargg = populate_Gamma(qqbargg); Omega_qqbargg = populate_OmegaI(qqbargg);
    Gamma_qgqg = populate_Gamma(qgqg); Omega_qgqg = populate_OmegaI(qgqg);
  }
  // charm quark tables
  if(hquark){
    Gamma_cqcq = populate_Gamma(cqcq); Omega_cqcq = populate_OmegaI(cqcq);
    Gamma_cgcg = populate_Gamma(cgcg); Omega_cgcg = populate_OmegaI(cgcg);
  }
}

scattering::~scattering()
{
}

// Determine if particles have engaged in a collision
//########################################################

bool scattering::isWounded(particle* p1, double dtau)
{
  int id = p1->getID();
  vector<double> fc_vel = p1->getFluidVelocity();
  vector<double> p1_lab = p1->getMomentum();
  vector<double> p1_cell = p1->lorentzboost(p1_lab,fc_vel);

  double E1_cell = p1_cell[0];
  double temp = p1->getFluidThermal()[2];
  double gamma = fc_vel[0]; 
  double dt = dtau/gamma;
  double prob = 1.0-exp(-incl_Gamma(E1_cell,temp,id)*dt/0.19733);
  double rand = drand48();
  
  if(rand < prob) return true;
  else return false;
} 

// Given that p1 scatters, determine collision partner
// p2 and collision products p3 and p4
//######################################################

vector<particle*> scattering::scatter(particle* p1)
{   
  // boost hard probe (p1) to fluid cell (fc) rest frame and rotate 

  vector<double> p1_x = p1->getPosition();
  vector<double> fc_vel = p1->getFluidVelocity();
  vector<double> p1_lab = p1->getMomentum();
  vector<double> p1_cell = p1->lorentzboost(p1_lab,fc_vel);
  double theta1 = p1->getTheta(p1_cell);
  double phi1 = p1->getPhi(p1_cell);
  p1_cell = p1->rotate(p1_cell,theta1,phi1); p1->setMomentum(p1_cell);
  if(abs(p1_cell[1]) > 0.000001 || abs(p1_cell[2]) > 0.000001) abort();
  
  // sample allowable scattering channels to determine ij -> kl, E2

  double temp = p1->getFluidThermal()[2];
  double db = get_cutoff(temp);
  double E1_cell = p1->getMomentum()[0];
  vector<int> ij2kl = sample_2to2(E1_cell,temp,p1->getID());

  // sample E2 and theta12 in the cell frame and fix kinematic s-limits

  particle* p2 = new particle(paraRdr); 
  p2->setID(ij2kl[1]); 
  p2->setPosition(p1_x);
  vector<double> p2_cell(4,0.0);
  double E2_cell = sample_E2(E1_cell,temp,ij2kl);
  double theta12;
  do{
    theta12 = sample_theta12(E1_cell,E2_cell,temp,ij2kl);
    if(theta12!=theta12){cout << "insufficient s" << endl; abort();}
  }while(theta12!=theta12);
  double phi12 = sample_phi12();
  p2_cell[0] = E2_cell;
  p2_cell[1] = E2_cell*sin(theta12)*cos(phi12);
  p2_cell[2] = E2_cell*sin(theta12)*sin(phi12);
  p2_cell[3] = E2_cell*cos(theta12);
  double slim0 = get_s_limits(E1_cell,E2_cell,db,ij2kl)[0];
  double sinit = get_MandelstamS(p1_cell,p2_cell);
  double tinit = get_MandelstamT(p1_cell,p2_cell); 
  double uinit = get_MandelstamU(p1_cell,p2_cell); 

  // boost p1 and p2 to their cms frame and rotate so p1 points towards zhat

  vector<double> u_cms = getUcms(p1_cell, p2_cell);
  vector<double> p1_cms = p1->lorentzboost(p1_cell, u_cms);
  vector<double> p2_cms = p2->lorentzboost(p2_cell, u_cms);
  double theta1_cms = p1->getTheta(p1_cms);
  double phi1_cms = p1->getPhi(p1_cms);
  p1_cms = p1->rotate(p1_cms,theta1_cms,phi1_cms);
  p2_cms = p2->rotate(p2_cms,theta1_cms,phi1_cms);
  p1->setMomentum(p1_cms);
  p2->setMomentum(p2_cms);

  // scatter p1 and p2, returning deflected theta13 and phi13

  double s = get_MandelstamS(p1_cms, p2_cms);
  double p1_mag = pow(p1_cms[1]*p1_cms[1]+p1_cms[2]*p1_cms[2]+p1_cms[3]*p1_cms[3],0.5);
  double theta13 = sample_theta13(s,p1_mag,temp,ij2kl);
  double phi13 = sample_phi13();

  // create deflected particles p3 and p4
  
  particle* p3 = new particle(paraRdr);
  vector<double> p3_cms;
  p3->setID(ij2kl[2]);
  p3->setPosition(p1_x);
  p3_cms = p3->inv_rotate(p1_cms, theta13, phi13);
  p3_cms = p3->inv_rotate(p3_cms,theta1_cms,phi1_cms);
  p3->setMomentum(p3_cms);

  particle* p4 = new particle(paraRdr);
  vector<double> p4_cms;
  p4->setID(ij2kl[3]);
  p4->setPosition(p1_x);
  p4_cms = p4->inv_rotate(p2_cms, theta13, phi13);
  p4_cms = p4->inv_rotate(p4_cms,theta1_cms,phi1_cms);
  p4->setMomentum(p4_cms);

  // restore p3 and p4 to the cell frame and then the lab frame

  vector<double> u_cms_inv = p3->reflectfourvector(u_cms); 
  vector<double> p3_cell = p3->lorentzboost(p3_cms, u_cms_inv);
  vector<double> p4_cell = p4->lorentzboost(p4_cms, u_cms_inv);
  p3_cell = p3->inv_rotate(p3_cell,theta1,phi1); p3->setMomentum(p3_cell); 
  p4_cell = p4->inv_rotate(p4_cell,theta1,phi1); p4->setMomentum(p4_cell);
  
  vector<double> u_cell_inv = p3->reflectfourvector(fc_vel);
  vector<double> p3_lab = p3->lorentzboost(p3_cell, u_cell_inv); p3->setMomentum(p3_lab);
  vector<double> p4_lab = p4->lorentzboost(p4_cell, u_cell_inv); p4->setMomentum(p4_lab);
   
  double sfinal = get_MandelstamS(p3_lab,p4_lab); 
  double tfinal = get_MandelstamT(p3_lab,p4_lab); 
  double ufinal = get_MandelstamU(p3_lab,p4_lab);
  double mandelstam_check = abs(sfinal-sinit) + abs(tfinal-tinit) + abs(ufinal-uinit);
  if(mandelstam_check > 0.000001){cout << "Mandelstam {s,t,u}: not conserved" << endl; abort();}

  // return the scattered products p3 and p4

  vector<particle*> products;
  products.push_back(p3); products.push_back(p4);
  delete p1, p2;
  return products;
}

// Sample the thermal target energy E2 (from dGamma/dE2)
//###########################################################################################

double scattering::sample_E2(double E1, double temp, vector<int> ij2kl)
{
  int npart = 500;
  double db = get_cutoff(temp);
  double dGdE2=0.0, dGdE2_max=0.0, dGdE2_rand;
  double E2 = E2_min, E2_rand;

  for(int n=0; n<npart; n++){
    dGdE2 = dGammadE2(E1,E2,temp,ij2kl);
    if(dGdE2 > dGdE2_max) dGdE2_max = dGdE2;
    E2 += dE2;
  }
  do{
    E2_rand = E2_min+drand48()*(E2_max-E2_min);
    dGdE2_rand = drand48()*dGdE2_max;
  }while(dGdE2_rand > dGammadE2(E1,E2_rand,temp,ij2kl));
  return E2_rand;
}

double scattering::dGammadE2(double E1, double E2, double temp, vector<int> ij2kl)
{
  double mass = 0.0;
  if(ij2kl[0]==4) mass = mc;
  double Q0 = get_cutoff(temp);
  vector<double> slim = get_s_limits(E1, E2, Q0, ij2kl);
  double s_min = slim[0], s_max = slim[1];
  double dGammadE2 = F(E2,temp,ij2kl[1])*Omega(E1,E2,Q0,ij2kl);
  return dGammadE2/(16.0*pow(M_PI,2.0)*E1*pow(E1*E1-mass*mass,0.5));
}

// Sample the thermal target angle theta12 (evaluated in the cell frame)
//#####################################################################################

double scattering::sample_theta12(double E1, double E2, double temp, vector<int> ij2kl)
{
  double mass = 0.0;
  if(ij2kl[0]==4) mass = mc;
  double Q0 = get_cutoff(temp);
  vector<double> slim = get_s_limits(E1, E2, Q0, ij2kl);
  double s_min = slim[0], s_max = slim[1], s, arg; 
  do{
    s = invert_OmegaI(E1,E2,temp,ij2kl,drand48());
    arg = (2.0*E1*E2+mass*mass-s)/(2.0*pow(E1*E1-mass*mass,0.5)*E2);
    if(arg!=arg || arg<-1 || arg>1) cout << "bad s" << endl;
  }while(arg!=arg || arg<-1 || arg>1);
  return acos(arg);
}

double scattering::sample_phi12()
{
  double phi = drand48()*2.0*M_PI;
  return phi;
}

// Sample angles theta & phi between p1 and p3 (formalism only valid in CMS frame)
//################################################################################

double scattering::sample_theta13(double s, double pmag, double temp, vector<int> ij2kl)
{
  double db = get_cutoff(temp);
  vector<double> tlim = get_t_limits(s, db, ij2kl);
  double t_min = tlim[0], t_max = tlim[1], t;
  do{
    t = invert_sigmaI(s,temp,ij2kl,drand48());
    if(t!=t || t<t_min || t>t_max) cout << "bad t" << endl;
  }while(t!=t || t<t_min || t>t_max);

  return acos(1.0 + t/(2.0*pmag*pmag)); 
}

double scattering::sample_phi13()
{
  double phi = drand48()*2.0*M_PI;
  return phi;
}

// Sample particle ID's {i,j,k,l} in 2 to 2 process ij -> kl
//#################################################################

vector<int> scattering::sample_2to2(double E1, double T, int i)
{
  vector<int> ij2kl(4,0);
  double rchan = drand48();
  int j, k, l;

  // gluon channel
 
  if(i==0){

    double gggg = interp_Gamma(E1,T,Gamma_gggg);
    double ggqqbar = interp_Gamma(E1,T,Gamma_ggqqbar);
    double gqgq = interp_Gamma(E1,T,Gamma_gqgq);
    double incl = gggg + ggqqbar + 6.0*gqgq;

    double p0 = 0.0;
    double p1 = gggg/incl;
    double p2 = (gggg+ggqqbar)/incl;
    double p3 = (gggg+ggqqbar+6.0*gqgq)/incl;

    if(abs(1.0-p3) > .00001){cout << "partition error: aborting." << endl; abort();}
    if(p0<=rchan && rchan<p1){
      j = 0;
      k = 0;
      l = 0;
    }
    else if(p1<=rchan && rchan<p2){
      j = 0;
      k = (rand()%3 + 1)*((rand()%2) * 2 - 1);
      l = -k;
    }
    else{
      j = (rand()%3 + 1)*((rand()%2) * 2 - 1);
      k = 0;
      l = j;
    }
  }
  // quark channel
 
  else if(abs(i)>0 && abs(i)<4){

    double qigqig = interp_Gamma(E1,T,Gamma_qgqg);
    double qiqjqiqj = interp_Gamma(E1,T,Gamma_qiqjqiqj);
    double qiqiqiqi = interp_Gamma(E1,T,Gamma_qiqiqiqi);
    double qiqibarqiqibar = interp_Gamma(E1,T,Gamma_qiqibarqiqibar);
    double qiqibarqjqjbar = interp_Gamma(E1,T,Gamma_qiqibarqjqjbar);
    double qiqibargg = interp_Gamma(E1,T,Gamma_qqbargg);
    double incl = qigqig + 4.0*qiqjqiqj + qiqiqiqi + qiqibarqiqibar + qiqibarqjqjbar + qiqibargg;

    double p0 = 0.0;
    double p1 = qigqig/incl;
    double p2 = (qigqig+4.0*qiqjqiqj)/incl;
    double p3 = (qigqig+4.0*qiqjqiqj+qiqiqiqi)/incl;
    double p4 = (qigqig+4.0*qiqjqiqj+qiqiqiqi+qiqibarqiqibar)/incl;
    double p5 = (qigqig+4.0*qiqjqiqj+qiqiqiqi+qiqibarqiqibar+qiqibarqjqjbar)/incl;
    double p6 = (qigqig+4.0*qiqjqiqj+qiqiqiqi+qiqibarqiqibar+qiqibarqjqjbar+qiqibargg)/incl;

    if(abs(1.0-p6) > .00001) {cout << "partition error: aborting." << endl; abort();}
    if(p0<=rchan && rchan<p1){
      j = 0;
      k = i;
      l = 0;
    }
    else if(p1<=rchan && rchan<p2){
      do{
	j =  (rand()%3 + 1)*((rand()%2) * 2 - 1);
      }while(abs(j)==abs(i));
      k = i;
      l = j;
    }
    else if(p2<=rchan && rchan<p3){
      j = i;
      k = i;
      l = i;
    }
    else if(p3<=rchan && rchan<p4){
      j = -i;
      k = i;
      l = -i;
    }
    else if(p4<=rchan && rchan<p5){
      j = -i;
      do{
	k =  (rand()%3 + 1)*((rand()%2) * 2 - 1);
      }while(abs(k)==abs(i));
      l = -k;
    }
    else{
      j = -i;
      k = 0;
      l = 0;
      }
  }

  // charm quark channel
 
  else if(i==4){

    double cqcq = 6.0*interp_Gamma(E1,T,Gamma_cqcq);
    double cgcg = interp_Gamma(E1,T,Gamma_cgcg);
    double incl = cqcq+cgcg;

    double p0 = 0.0;
    double p1 = cqcq/incl;
    double p2 = (cqcq+cgcg)/incl;
    if(abs(1.0-p2) > .00001){cout << "partition error: aborting." << endl; abort();}

    if(p0<=rchan && rchan<p1){
      j = (int(rand()%3) + 1)*(int(rand()%2)*2 - 1);
      k = 4;
      l = j;
    }
    else{
      j = 0;
      k = 4;
      l = 0;
    }
  }
  else{
    cout << "no such particle i" << endl;
  }

  ij2kl[0] = i; ij2kl[1] = j; ij2kl[2] = k; ij2kl[3] = l;
  return ij2kl;
}

// Determine inclusive scattering rate for particle species "i" 
//############################################################

double scattering::incl_Gamma(double E1, double T, int i)
{
  if(i==0){ 
    return interp_Gamma(E1,T,Gamma_gggg) + interp_Gamma(E1,T,Gamma_ggqqbar) + 6.0*interp_Gamma(E1,T,Gamma_gqgq);
  }
  else if(abs(i)>0 && abs(i)<4){
    return interp_Gamma(E1,T,Gamma_qgqg) + 4.0*interp_Gamma(E1,T,Gamma_qiqjqiqj) + interp_Gamma(E1,T,Gamma_qiqiqiqi) + 
           interp_Gamma(E1,T,Gamma_qiqibarqiqibar) + interp_Gamma(E1,T,Gamma_qiqibarqjqjbar) + interp_Gamma(E1,T,Gamma_qqbargg);
  }
  else if(i==4){
    return 6.0*interp_Gamma(E1,T,Gamma_cqcq) + interp_Gamma(E1,T,Gamma_cgcg);
  }
  else{
    cout << "warning: disallowed particle 'i' in ij->kl channel" << endl;
  }
}

// Determine scattering rate for channel ij->kl
//##############################################################

double scattering::Gamma(double E1, double T, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

  // light processes
  
  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return interp_Gamma(E1,T,Gamma_gggg);
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    return interp_Gamma(E1,T,Gamma_ggqqbar);
  }
  // {g q -> g q, g qbar -> g qbar}
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    return interp_Gamma(E1,T,Gamma_qgqg);
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){
    return interp_Gamma(E1,T,Gamma_qiqjqiqj);
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){
    return interp_Gamma(E1,T,Gamma_qiqiqiqi);
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    return interp_Gamma(E1,T,Gamma_qiqibarqjqjbar);
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    return interp_Gamma(E1,T,Gamma_qiqibarqiqibar);
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    return interp_Gamma(E1,T,Gamma_qqbargg);
  }

  // heavy processes
  
  // {cq -> cq}
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    return interp_Gamma(E1,T,Gamma_cqcq);
  }
  // {cg -> cg}
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    return interp_Gamma(E1,T,Gamma_cgcg);
  }

  // >>no such event<<

  else return 0.0;
}

// Interpolate scattering rate table for arbitary E1 and T
//##################################################################

double scattering::interp_Gamma(double E1, double T, double **Gamma)
{
  // see Wikipedia bilinear interpolation

  int iE1 = (E1-E1_min)/double(dE1);
  int iE1p1;
  if(nE1>1) iE1p1 = (E1-E1_min)/dE1+1;
  else iE1p1 = (E1-E1_min)/dE1;
  double x1 = E1_min + iE1*dE1;
  double x2 = E1_min + (iE1+1)*dE1;
  double x = E1;

  int iT = (T-T_min)/dT;
  int iTp1;
  if(nT>1) iTp1 = (T-T_min)/dT+1;
  else iTp1 = (T-T_min)/dT;
  double y1 = T_min + iT*dT;
  double y2 = T_min + (iT+1)*dT;
  double y = T;
  double F11 = Gamma[iE1][iT];
  double F12 = Gamma[iE1][iTp1];
  double F21 = Gamma[iE1p1][iT];
  double F22 = Gamma[iE1p1][iTp1];

  if(E1_min==E1_max && T_min==T_max) return F11;
  else if(E1_min==E1_max) return (F11*y - F12*y + F12*y1 - F11*y2)/(y1 - y2);
  else if(T_min==T_max) return (F11*x - F21*x + F21*x1 - F11*x2)/(x1 - x2);
  else return (F11*(x2-x)*(y2-y)+F21*(x-x1)*(y2-y)+F12*(x2-x)*(y-y1)+F22*(x-x1)*(y-y1))/((x2-x1)*(y2-y1));
}

// Determine scattering rate tables for all channels ij->kl
//#####################################################################################

double** scattering::populate_Gamma(vector<int> ij2kl)
{
  double** Gamma;
  int counter=0;

  Gamma = new double*[nE1];
  for(int iE1=0; iE1<nE1; iE1++){
    double E1 = E1_min + iE1*dE1;
    Gamma[iE1] = new double[nT]();
    for(int iT=0; iT<nT; iT++){ 
      counter++;
      double T = T_min + iT*dT;
      Gamma[iE1][iT] = calc_Gamma(E1,T,ij2kl);
      drawProgressBar(60,double(counter)/double(nE1*nT));
    }
  }

  return Gamma;
}

// Integration to calculate the scattering rate for a specific channel ij->kl
//############################################################################

double scattering::calc_Gamma(double E1, double T, vector<int> ij2kl)
{  
  // integration precision
  int ns = 500;
  
  double Gamma = 0.0, sigma=0.0, Omega=0.0;
  double mass = 0.0;
  if(ij2kl[0]==4) mass = mc;
  double db = get_cutoff(T);
      
  for(int iE2=0; iE2<nE2; iE2++){// E2
    double E2 = E2_min + iE2*dE2;
    vector<double> slim = get_s_limits(E1,E2,db,ij2kl);
    double s_min = slim[0], s_max = slim[1], ds = (s_max-s_min)/double(ns);

    for(int is=0; is<ns; is++){ // s 
      double s = s_min + is*ds;
      vector<double> tlim = get_t_limits(s,db,ij2kl); 
      double t_min = tlim[0], t_max = tlim[1];
      sigma = max(sigmaI(s,t_max,T,ij2kl)-sigmaI(s,t_min,T,ij2kl),0.0);
      Omega += (s-mass*mass)*sigma*ds;
      Gamma += F(E2,T,ij2kl[1])*(s-mass*mass)*sigma/(16.0*pow(M_PI,2.0)*E1*pow(E1*E1-mass*mass,0.5))*ds*dE2;
    } // s
  } // E2

  return Gamma;
}


// Integral method to invert dsigma / dt and sample t (for theta13)
// ######################################################################################

double scattering::invert_sigmaI(double s, double temp, vector<int> ij2kl, double frac)
{
  int counts=0;
  double db = get_cutoff(temp);
  vector<double> tlim = get_t_limits(s, db, ij2kl);
  double t_min = tlim[0], t_max = tlim[1];

  double t_lb = t_min, t_ub = t_max, t = t_min+drand48()*(t_max-t_min);
  double sigma = (sigmaI(s,t_max,temp,ij2kl)-sigmaI(s,t_min,temp,ij2kl));
  double sigma_target = sigmaI(s,t_min,temp,ij2kl)+frac*sigma;  
  double dist = sigma, precision = 0.001*sigma;

  do{
    dist = (sigmaI(s,t,temp,ij2kl) - sigma_target);
    if(dist<0.0){
      t_lb = t;
      t = t+drand48()*(t_ub-t);
    }
    else{
      t_ub = t;
      t = t_lb+drand48()*(t-t_lb);
    }
    counts++;
    if(counts>50){cout << "invert_sigma: exceeded maximum repetitions. check precision." << endl; abort();} 
  }while(abs(dist) > precision);

  return t;
}

// Used to sample dsigma / dt and sample t (for theta13)
// Indefinite sigma_ij->kl(s,t) = 1/(1+delta_kl)  1/(16 pi (s-mc^2)^2)  \int  dt Msq_ij->kl(s)
//#######################################################################################

double scattering::sigmaI(double s, double t, double temp, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];
  double db = get_cutoff(temp);

  // light processes

  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return (9.0*pow(gs,4.0)*(-(pow(s,4.0)/t) + 3.0*pow(s,2.0)*t + (s*pow(t,2.0))/2.0 + pow(t,3.0)/3.0 - pow(s,4.0)/(s+t) + pow(s,3.0)*log(abs(t)) - pow(s,3.0)*log(abs(s+t))))/(64.0*M_PI*pow(s,4.0));
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    return -((pow(gs,4.0)*(17.0*pow(s,2.0)*t + 9.0*s*pow(t,2.0) + 6.0*pow(t,3.0) + 4.0*pow(s,3.0)*log(abs(t)) - 4.0*pow(s,3.0)*log(abs(s+t))))/(384.0*M_PI*pow(s,4.0)));
  }
  // {g q -> g q, g qbar -> g qbar}
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    return (pow(gs,4.0)*(-((18.0*pow(s,3.0))/t) + 13.0*s*t + 2.0*pow(t,2.0) + 18.0*pow(s,2.0)*log(abs(t)) + 4.0*pow(s,2.0)*log(abs(s+t))))/(144.0*M_PI*pow(s,3.0));
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){
    return 4.0/9.0*pow(gs,4.0)*(-((2.0*pow(s,2.0))/t) + t + 2*s*log(abs(t)))/(16.0*M_PI*pow(s,2.0));
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){
    return 8.0/27.0*pow(gs,4.0)*(-((3.0*pow(s,2.0))/t) + 3.0*t - (3.0*pow(s,2.0))/(s+t) + 4.0*s*log(abs(t)) - 4.0*s*log(abs(s+t)))/(32.0*M_PI*pow(s,2.0));
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    return (4.0*pow(gs,4.0)*(pow(s,2.0)*t + s*pow(t,2.0) + (2.0*pow(t,3.0))/3.0))/(144.0*M_PI*pow(s,4.0));
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    return (8.0*pow(gs,4.0)*(-((3.0*pow(s,4.0))/t) + pow(s,2.0)*t + s*pow(t,2.0) + pow(t,3.0) + 2.0*pow(s,3.0)*log(abs(t))))/(432.0*M_PI*pow(s,4.0));
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    return -((8.0*pow(gs,4.0)*(17.0*pow(s,2.0)*t + 9.0*s*pow(t,2.0) + 6.0*pow(t,3.0) + 4.0*pow(s,3.0)*log(abs(t)) - 4.0*pow(s,3.0)*log(abs(s+t))))/(864.0*M_PI*pow(s,4.0)));
  }

  // heavy processes
  
  // {cq -> cq}
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    return (pow(gs,4)*(pow(db,4) + 10*pow(mc,4) - 8*pow(mc,2)*s + 2*pow(s,2) - pow(t,2) + pow(db,2)*(-4*pow(mc,2) + 2*s + t) - 2*(pow(db,2) + 2*pow(mc,2) - s)*(pow(db,2) - t)*log(pow(db,2) - t)))/
           (36.*M_PI*pow(pow(mc,2) - s,2)*(pow(db,2) - t));
  }
  // {cg -> cg}
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    return (pow(gs,4)*((-18*pow(pow(mc,2) - s,3)*(pow(db,2) - pow(mc,2) + s))/(pow(db,2) - t) + (13*pow(mc,4) - 10*pow(mc,2)*s + 13*pow(s,2))*t + 2*(-pow(mc,2) + s)*pow(t,2) - 
           (16*pow(mc,4)*pow(pow(mc,2) - s,2))/(-pow(mc,2) + s + t) - (18*pow(mc,2)*pow(pow(mc,2) - s,3)*log(-pow(mc,2) + s + t))/(db - pow(mc,2) + s) + ((pow(mc,2) - s)*(9*(pow(db,2) - 
           2*pow(mc,2) + 2*s)*(-pow(mc,4) + s*(pow(db,2) + s))*log(pow(db,2) - t) + (-42*pow(mc,6) + 84*pow(mc,4)*s - 46*pow(mc,2)*pow(s,2) + 4*pow(s,3) + pow(db,2)*(15*pow(mc,4) - 
	   15*pow(mc,2)*s + 4*pow(s,2)))*log(-pow(mc,2) + s + t)))/(-pow(db,2) + pow(mc,2) - s)))/(144.*M_PI*pow(pow(mc,2) - s,4));
  }

  // >>no such event<<

  else{
    return 0.0;
  }
}

// Calculate Omega from the Omega interpolation table
//##########################################################################################

double scattering::Omega(double E1, double E2, double temp, vector<int> ij2kl)
{
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];
  double db = get_cutoff(temp); 
  vector<double> slim = get_s_limits(E1,E2,db,ij2kl);
  double s_min = slim[0], s_max = slim[1];

  // light processes
  
  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return interp_OmegaI(s_max,temp,ij2kl,Omega_gggg)-interp_OmegaI(s_min,temp,ij2kl,Omega_gggg);
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_ggqqbar)-interp_OmegaI(s_min,temp,ij2kl,Omega_ggqqbar);
  }
  // {g q -> g q, g qbar -> g qbar}
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qgqg)-interp_OmegaI(s_min,temp,ij2kl,Omega_qgqg);
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qiqjqiqj)-interp_OmegaI(s_min,temp,ij2kl,Omega_qiqjqiqj);
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qiqiqiqi)-interp_OmegaI(s_min,temp,ij2kl,Omega_qiqiqiqi);
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qiqibarqjqjbar)-interp_OmegaI(s_min,temp,ij2kl,Omega_qiqibarqjqjbar);
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qiqibarqiqibar)-interp_OmegaI(s_min,temp,ij2kl,Omega_qiqibarqiqibar);
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_qqbargg)-interp_OmegaI(s_min,temp,ij2kl,Omega_qqbargg);
  }

  // heavy processes
  
  // {cq -> cq}
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_cqcq)-interp_OmegaI(s_min,temp,ij2kl,Omega_cqcq);
  }
  // {cg -> cg}
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    return interp_OmegaI(s_max,temp,ij2kl,Omega_cgcg)-interp_OmegaI(s_min,temp,ij2kl,Omega_cgcg);
  }

  // >>no such event<<

  else return 0.0;
}

// Interpolate the OmegaI table for arbitrary Mandelstam-s
//###################################################################

double scattering::interp_OmegaI(double s, double temp, vector<int> ij2kl, double **OmegaI)
{
  // see Wikipedia bilinear interpolation

  // table integration precision
  int ns = ns_table;
  double mass = 0.0; 
  if(ij2kl[0]==4) mass = mc;
  double s_min = mass*mass, s_max = pow(E1_max+E2_max,2.0), ds = (s_max-s_min)/double(ns);

  int is = (s-s_min)/double(ds);
  int isp1;
  if(ns>1) isp1 = (s-s_min)/ds+1;
  else isp1 = (s-s_min)/ds;
  double x1 = s_min + is*ds;
  double x2 = s_min + (is+1)*ds;
  double x = s;
  
  double T = temp;
  int iT = (T-T_min)/dT;
  int iTp1;
  if(nT>1) iTp1 = (T-T_min)/dT+1;
  else iTp1 = (T-T_min)/dT;
  double y1 = T_min + iT*dT;
  double y2 = T_min + (iT+1)*dT;
  double y = T;

  double F11 = OmegaI[is][iT];
  double F12 = OmegaI[is][iTp1];
  double F21 = OmegaI[isp1][iT];
  double F22 = OmegaI[isp1][iTp1];

  if(s_min==s_max && T_min==T_max) return F11;
  else if(s_min==s_max) return (F11*y - F12*y + F12*y1 - F11*y2)/(y1 - y2);
  else if(T_min==T_max) return (F11*x - F21*x + F21*x1 - F11*x2)/(x1 - x2);
  else return (F11*(x2-x)*(y2-y)+F21*(x-x1)*(y2-y)+F12*(x2-x)*(y-y1)+F22*(x-x1)*(y-y1))/((x2-x1)*(y2-y1));
}

// Populate the Omega table
//####################################################################

double** scattering::populate_OmegaI(vector<int> ij2kl)
{
  double** OmegaI;
  int counter=0;
  int ns = ns_table;
  double mass = 0.0; 
  if(ij2kl[0]==4) mass = mc;
  double s_min = mass*mass, s_max = pow(E1_max+E2_max,2.0), ds = (s_max-s_min)/double(ns);

  OmegaI = new double*[ns];
  for(int is=0; is<ns; is++){ 
    double s = s_min + is*ds;
    OmegaI[is] = new double[nT]();
    for(int iT=0; iT<nT; iT++){ 
      double T = T_min + iT*dT;
      OmegaI[is][iT] = calc_OmegaI(s,T,ij2kl);
      drawProgressBar(60,double(counter)/double(ns*nT)); 
      counter++;
    }
  }

  return OmegaI;
}

// Integrate (s-mc^2) sigma(s) from 0 to s to calc OmegaI
//####################################################################

double scattering::calc_OmegaI(double s, double temp, vector<int> ij2kl)
{
  // integration precision
  int ns = 200;

  double mass = 0.0; 
  if(ij2kl[0]==4) mass = mc;
  double s_min = mass*mass, s_max = s, ds = (s_max-s_min)/double(ns);
  double db = get_cutoff(temp);
  double sigma = 0.0, OmegaI = 0.0;

  for(int is=0; is<ns; is++){ // s 
    double ss = s_min+is*ds;
    vector<double> tlim = get_t_limits(ss,db,ij2kl); 
    double t_min = tlim[0], t_max = tlim[1];
    sigma = max(sigmaI(ss,t_max,temp,ij2kl)-sigmaI(ss,t_min,temp,ij2kl),0.0);
    if(sigma!=sigma) sigma = 0.0;
    OmegaI += (ss-mass*mass)*sigma*ds;
  } // s

  return OmegaI;
}

// Invert \int{0,s} (s-m*m)*sigma(s) to sample s
//################################################################################################## 

double scattering::invert_OmegaI(double E1, double E2, double temp, vector<int> ij2kl, double frac)
{
  int counts=0;
  double db = get_cutoff(temp);
  vector<double> slim = get_s_limits(E1, E2, db, ij2kl);
  double s_min = slim[0], s_max = slim[1];
  double s_lb = s_min, s_ub = s_max, s = s_min+drand48()*(s_max-s_min);
  double Omega0 = (OmegaI(s_max,temp,ij2kl)-OmegaI(s_min,temp,ij2kl));
  double Omega_target = OmegaI(s_min,temp,ij2kl)+frac*Omega0; 
  double dist = Omega0, precision = 0.001*Omega0;
  
  do{
    dist = (OmegaI(s,temp,ij2kl) - Omega_target);
    if(dist<0.0){
      s_lb = s;
      s = (s+s_ub)/2.0;
    }
    else{
      s_ub = s;
      s = (s_lb+s)/2.0;
    }
    counts++;
    if(counts>50){cout << "invert_Omega: exceeded maximum repetitions. check precision." << endl; abort();} 
  }while(abs(dist) > precision); 
  return s;
}

// indefinite Omega_ij->kl(s, Q0) = \int{0,s} (s-m*m)*sigma(s)
//############################################################################################

double scattering::OmegaI(double s, double temp, vector<int> ij2kl)
{
  double Q0 = get_cutoff(temp);
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];
  
  // light processes

  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){ 
    return interp_OmegaI(s, temp, ij2kl, Omega_gggg);
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    return interp_OmegaI(s, temp, ij2kl, Omega_ggqqbar);
  }
  // {g q -> g q, g qbar -> g qbar}
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    return interp_OmegaI(s, temp, ij2kl, Omega_gqgq);
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){  
    return interp_OmegaI(s, temp, ij2kl, Omega_qiqjqiqj);
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){         
    return interp_OmegaI(s, temp, ij2kl, Omega_qiqiqiqi);
  }
  // {qi qibar -> qj qjbar} 
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    return interp_OmegaI(s, temp, ij2kl, Omega_qiqibarqjqjbar);
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    return interp_OmegaI(s, temp, ij2kl, Omega_cqcq);
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    return interp_OmegaI(s, temp, ij2kl, Omega_qqbargg);
  }

  // heavy processes

  // {cq -> cq}
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    return interp_OmegaI(s, temp, ij2kl, Omega_cqcq);
  }
  // {cg -> cg} 
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    return interp_OmegaI(s, temp, ij2kl, Omega_cgcg);
  }
  // >>no such event<<
  else{
    return 0.0;
  }
}

// Set the Debye screen mass
//###################################################################

double scattering::get_cutoff(double temp)
{
  return sqrt(3.0/2.0)*gs*temp;
}

// Min/max limits on Mandelstam-s
//########################################################################################

vector<double> scattering::get_s_limits(double E1, double E2, double db, vector<int> ij2kl)
{ 
  // **Evaluated for E1 and E2 in the fluid cell rest frame**
  vector<double> slim(2,0.0);
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){
    slim[0] = 2.0*db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {g g -> q qbar}
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    slim[0] = 2.0*db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {g q -> g q, g qbar -> g qbar}
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    slim[0] = 2.0*db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){
    slim[0] = db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){
    slim[0] = 2.0*db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    slim[0] = db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    slim[0] = db*db;
    slim[1] = 4.0*E1*E2;
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    slim[0] = 2.0*db*db;
    slim[1] = 4.0*E1*E2;
  }
  // cg -> cg 
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    slim[0] = mc*mc + 2.0*E1*E2 - 2.0*sqrt(E1*E1-mc*mc)*E2;
    slim[1] = mc*mc + 2.0*E1*E2 + 2.0*sqrt(E1*E1-mc*mc)*E2;
  }
  // cq -> cq 
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    slim[0] = mc*mc + 2.0*E1*E2 - 2.0*sqrt(E1*E1-mc*mc)*E2;
    slim[1] = mc*mc + 2.0*E1*E2 + 2.0*sqrt(E1*E1-mc*mc)*E2;
  }
  else{
    cout << "no such process" << endl;
    abort();
  }
  return slim;
}

// Min/max limits on Mandelstam-t
//########################################################################################


vector<double> scattering::get_t_limits(double s, double db, vector<int> ij2kl)
{
  vector<double> tlim(2,0.0);
  int i = ij2kl[0], j = ij2kl[1], k = ij2kl[2], l = ij2kl[3];

  // {g g -> g g} 
  if(i==0 && j==0 && k==0 && l==0){
    tlim[0] = -s+db*db;
    tlim[1] = -db*db;
  }
  // {g g -> q qbar} 
  else if(i==0 && j==0 && abs(k)>0 && abs(k)<4 && l==-k){
    tlim[0] = -s+db*db;
    tlim[1] = -db*db;
  }
  // {g q -> g q, g qbar -> g qbar} 
  else if((i==0 && abs(j)>0 && abs(j)<4 && k==0 && l==j) || (abs(i)>0 && abs(i)<4 && j==0 && k==i && l==0)){
    tlim[0] = -s+db*db;
    tlim[1] = -db*db;
  }
  // {qi qj -> qi qj, qi qjbar -> qi qjbar, qibar qj -> qibar qj, qibar qjbar -> qibar qjbar} 
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==k && j==l && i!=j && k!=l){
    tlim[0] = -s;
    tlim[1] = -db*db;
  }
  // {qi qi -> qi qi, qibar qibar -> qibar qibar} 
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==j && j==k && k==l){
    tlim[0] = -s+db*db;
    tlim[1] = -db*db;
  }
  // {qi qibar -> qj qjbar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i!=k && j!=l){
    tlim[0] = -s;
    tlim[1] = -db*db;
  }
  // {qi qibar -> qi qibar}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && abs(k)>0 && abs(k)<4 && abs(l)>0 && abs(l)<4 && i==-j && k==-l && i==k && j==l){
    tlim[0] = -s;
    tlim[1] = -db*db;
  }
  // {q qbar -> gg}
  else if(abs(i)>0 && abs(i)<4 && abs(j)>0 && abs(j)<4 && k==0 && l==0 && i==-j){
    tlim[0] = -s+db*db;
    tlim[1] = -db*db;
  }
  // cg -> cg 
  else if((i==0 && j==4 && k==0 && l==4) || (i==4 && j==0 && k==4 && l==0)){
    tlim[0] = -pow(s-mc*mc,2.0)/s; 
    tlim[1] = 0.0;
  }
  // cq -> cq 
  else if((abs(i)>0 && abs(i)<4 && abs(j)==4 && k==i && l==j) || (abs(i)==4 && abs(j)>0 && abs(j)<4 && k==i && l==j)){
    tlim[0] = -pow(s-mc*mc,2.0)/s; 
    tlim[1] = 0.0;
  }
  else{
    cout << "no such channel" << endl;
    abort();
  }
  return tlim;
}

// Fermi-Dirac/Bose-Einstein distribution function 
//####################################################################

double scattering::F(double E2, double temp, int j)
{
  if(j==0) return gg*pow(exp(E2/temp)-1.0,-1.0);
  else return gq*pow(exp(E2/temp)+1.0,-1.0);
}

// Get the Lorentz invariant Mandelstam variables
//#######################################################################

double scattering::get_MandelstamS(vector<double> &p1, vector<double> &p2)
{
  double v0, v1, v2, v3, s;
  v0 = pow(p1[0]+p2[0],2.0);
  v1 = pow(p1[1]+p2[1],2.0);
  v2 = pow(p1[2]+p2[2],2.0);
  v3 = pow(p1[3]+p2[3],2.0);
  s = v0-v1-v2-v3;
  return s;
}

double scattering::get_MandelstamT(vector<double> &p1, vector<double> &p3)
{
  double v0, v1, v2, v3, t;
  v0 = pow(p1[0]-p3[0],2.0);
  v1 = pow(p1[1]-p3[1],2.0);
  v2 = pow(p1[2]-p3[2],2.0);
  v3 = pow(p1[3]-p3[3],2.0);
  t = v0-v1-v2-v3;
  return t;
}

double scattering::get_MandelstamU(vector<double> &p1, vector<double> &p4)
{
  double v0, v1, v2, v3, u;
  v0 = pow(p1[0]-p4[0],2.0);
  v1 = pow(p1[1]-p4[1],2.0);
  v2 = pow(p1[2]-p4[2],2.0);
  v3 = pow(p1[3]-p4[3],2.0);
  u = v0-v1-v2-v3;
  return u;
}
