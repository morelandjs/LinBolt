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
  static vector<int> dim; 
  static vector<double> gridmin; 
  static vector<double> gridmax;
 public:
  static vector<int> getDim(){return dim;}
  static vector<double> getGridMin(){return gridmin;}
  static vector<double> getGridMax(){return gridmax;}
  void setDim(vector<int> _dim){dim = _dim;}
  void setGridMin(vector<double> _gridmin){gridmin = _gridmin;}
  void setGridMax(vector<double> _gridmax){gridmax = _gridmax;}
  vector<double> coord2position(vector<int> &_coord);
  vector<int> position2coord(vector<double> &_position);
  vector<double> lorentzboost(vector<double> &_fourvector, vector<double> &_boost);
  template<class type>
    type fourproduct(vector<type> v1, vector<type> v2){
    double product=0.0;
    product += v1[0]*v2[0];
    for(int i=1; i<4; i++){product -= v1[i]*v2[i];}
    return product;
  }
};

#endif
