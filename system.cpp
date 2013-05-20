#include "system.h"
#include "arsenal.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
using std::vector;
using namespace boost::assign;
using namespace std;

vector<int> system::dim;
vector<double> system::gridmin;
vector<double> system::gridmax;

vector<double> system::coord2position(vector<int> &_grid_coord)
{
  vector<double> position;
  for(int i=0; i<dim.size(); i++){
    position.push_back(gridmin[i]+_grid_coord[i]*(gridmax[i]-gridmin[i])/(dim[i]-1));
  }
  return position;
}

vector<int> system::position2coord(vector<double> &_four_position)
{
  vector<int> coordinates;
  for(int i=0; i<dim.size(); i++){
    coordinates.push_back(int((_four_position[i]-gridmin[i])/(gridmax[i]-gridmin[i])*(dim[i]-1)+0.5));
  }
  return coordinates;
}

vector<double> system::lorentzboost(vector<double> &_fourvector, vector<double> &_boost)
{
  vector<double> u = _fourvector;
  vector<double> v = _boost;
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

  vector<double> uprime;  
  uprime.push_back(gamma*u[0]-gamma*beta_x*u[1]-gamma*beta_y*u[2]-gamma*beta_z*u[3]);
  uprime.push_back(-gamma*beta_x*u[0]+(1+(gamma-1)*beta_xx)*u[1]+((gamma-1)*beta_xy)*u[2]+((gamma-1)*beta_xz)*u[3]);
  uprime.push_back(-gamma*beta_y*u[0]+((gamma-1)*beta_yx)*u[1]+(1+(gamma-1)*beta_yy)*u[2]+((gamma-1)*beta_yz)*u[3]);
  uprime.push_back(-gamma*beta_z*u[0]+((gamma-1)*beta_zx)*u[1]+((gamma-1)*beta_zy)*u[2]+(1+(gamma-1)*beta_zz)*u[3]);
  
  return uprime;
}
