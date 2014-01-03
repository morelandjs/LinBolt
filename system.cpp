#include "system.h"
#include "arsenal.h"
#include <iostream>
#include <cmath>
#include <vector>
using std::vector;
using namespace std;

vector<int> system::steps;
vector<double> system::gridmin;
vector<double> system::gridmax;

vector<double> system::coord2position(vector<int> &coord)
{
  vector<double> position(coord.size(),0.0);
  for(int i=0; i<coord.size(); i++){
    double ds = (gridmax[i]-gridmin[i])/max(steps[i]-1,1);
    position[i] = gridmin[i]+coord[i]*ds;
  }
  return position;
}

vector<int> system::position2coord(vector<double> &position)
{
  vector<int> coord(position.size(),0);
  for(int i=0; i<position.size(); i++){
    double frac = (position[i]-gridmin[i])/max(gridmax[i]-gridmin[i],0.0001);
    coord[i] = int(frac*steps[i]);
  }
  return coord;
}

double system::getTheta(vector<double> &v)
{
  double z = v[3];
  double r = pow(v[1]*v[1]+v[2]*v[2]+v[3]*v[3],0.5);
  return acos(z/r);
}

double system::getPhi(vector<double> &v)
{
  double x = v[1], y = v[2];
  double phi = atan2(y,x);
  return phi;
}

vector<double> system::rotate(vector<double> &v, double theta, double phi)
{
  vector<double> rv(4,0.0);
  rv[0] = 1.0*v[0];
  rv[1] = cos(theta)*cos(phi)*v[1]+cos(theta)*sin(phi)*v[2]-sin(theta)*v[3];
  rv[2] = -sin(phi)*v[1]+cos(phi)*v[2];
  rv[3] = sin(theta)*cos(phi)*v[1]+sin(theta)*sin(phi)*v[2]+cos(theta)*v[3];
  return rv;
}

vector<double> system::inv_rotate(vector<double> &v, double theta, double phi)
{
  vector<double> rv(4,0.0);
  rv[0] = 1.0*v[0];
  rv[1] = cos(theta)*cos(phi)*v[1]-sin(phi)*v[2]+sin(theta)*cos(phi)*v[3];
  rv[2] = cos(theta)*sin(phi)*v[1]+cos(phi)*v[2]+sin(theta)*sin(phi)*v[3];
  rv[3] = -sin(theta)*v[1]+cos(theta)*v[3];
  return rv;
}

vector<double> system::lorentzboost(vector<double> &u, vector<double> &_v)
{
  vector<double> v = _v;
  double gamma = v[0];
  double beta = pow(1-1/(gamma*gamma),0.5);
  double beta_x = v[1]/gamma;
  double beta_y = v[2]/gamma;
  double beta_z = v[3]/gamma;
  double beta_xx, beta_yy, beta_zz, beta_xy, beta_xz, beta_yx, beta_yz, beta_zx, beta_zy;
  if(beta<0.0000000001){beta_xx=1; beta_yy=1; beta_zz=1; beta_xy=1; beta_xz=1; beta_yx=1; beta_yz=1; beta_zx=1; beta_zy=1;}
  else{ 
    beta_xx = (beta_x*beta_x)/(beta*beta);
    beta_yy = (beta_y*beta_y)/(beta*beta);
    beta_zz = (beta_z*beta_z)/(beta*beta);
    beta_xy = (beta_x*beta_y)/(beta*beta);
    beta_xz = (beta_x*beta_z)/(beta*beta);
    beta_yx = (beta_y*beta_x)/(beta*beta);
    beta_yz = (beta_y*beta_z)/(beta*beta);
    beta_zx = (beta_z*beta_x)/(beta*beta);
    beta_zy = (beta_z*beta_y)/(beta*beta);
  }
  vector<double> uprime(4,0.0);  
  uprime[0] = (gamma*u[0]-gamma*beta_x*u[1]-gamma*beta_y*u[2]-gamma*beta_z*u[3]);
  uprime[1] = (-gamma*beta_x*u[0]+(1+(gamma-1)*beta_xx)*u[1]+((gamma-1)*beta_xy)*u[2]+((gamma-1)*beta_xz)*u[3]);
  uprime[2] = (-gamma*beta_y*u[0]+((gamma-1)*beta_yx)*u[1]+(1+(gamma-1)*beta_yy)*u[2]+((gamma-1)*beta_yz)*u[3]);
  uprime[3] = (-gamma*beta_z*u[0]+((gamma-1)*beta_zx)*u[1]+((gamma-1)*beta_zy)*u[2]+(1+(gamma-1)*beta_zz)*u[3]);
  return uprime;
}

vector<double> system::getUcms(vector<double> &p1, vector<double> &p2)
{
  double beta, betax, betay, betaz, gamma;
  vector<double> u_cms(4,0.0);
  betax = (p1[1]+p2[1])/(p1[0]+p2[0]);
  betay = (p1[2]+p2[2])/(p1[0]+p2[0]);
  betaz = (p1[3]+p2[3])/(p1[0]+p2[0]);
  beta = pow(betax*betax+betay*betay+betaz*betaz,0.5);
  gamma = 1/pow(1-beta*beta,0.5);
  u_cms[0] = gamma;
  u_cms[1] = gamma*betax;
  u_cms[2] = gamma*betay;
  u_cms[3] = gamma*betaz;
  return u_cms;
}
