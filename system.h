#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "ParameterReader.h"
#define PI 3.14159
#define clight 1.0
using std::vector;

class system
{
 protected:
  static vector<int> steps; 
  static vector<double> gridmin; 
  static vector<double> gridmax;
 public:
  static vector<int> getSteps(){return steps;}
  static vector<double> getGridMin(){return gridmin;}
  static vector<double> getGridMax(){return gridmax;}
  void setSteps(vector<int> _steps){steps = _steps;}
  void setGridMin(vector<double> _gridmin){gridmin = _gridmin;}
  void setGridMax(vector<double> _gridmax){gridmax = _gridmax;}
  vector<double> coord2position(vector<int> &coord);
  vector<int> position2coord(vector<double> &position);

  double getTheta(vector<double> &v);
  double getPhi(vector<double> &v);
  vector<double> rotate(vector<double> &v, double theta, double phi);
  vector<double> inv_rotate(vector<double> &v, double theta, double phi);
  vector<double> lorentzboost(vector<double> &u, vector<double> &v);
  vector<double> getUcms(vector<double> &p1, vector<double> &p2);

  template<class type>
    vector<type> reflectfourvector(vector<type> &v){
    vector<type> rv(4,0.0);
    rv[0] = v[0];
    for(int i=1; i<4; i++){rv[i] = -1.0*v[i];}
    return rv;
  }

  template<class type>
    type fourproduct(vector<type> &v1, vector<type> &v2){
    double product=0.0;
    product += v1[0]*v2[0];
    for(int i=1; i<4; i++){product -= v1[i]*v2[i];}
    return product;
  }
};

#endif
