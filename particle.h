#ifndef PARTICLE_H
#define PARTICLE_H

#include "system.h"
#include "arsenal.h"
#include "medium.h"
#include "ParameterReader.h"
#include <vector>
#include <map>
using std::vector;
using std::map;
using namespace std;

class particle : public system
{

 private:
  ParameterReader* paraRdr;

  map<int,double> masstable;
  double alphas, gs, pathlength, Qperp;

  int id;
  vector<int> coord;
  vector<double> position;
  vector<double> velocity;
  vector<double> momentum;
  vector<double> fluid_velocity;
  vector<double> fluid_thermal;

 public:
  // constructor/de-constructor
  particle(ParameterReader* _paraRdr);
  ~particle();

  // get functions
  int getID(){return id;}
  vector<int> getCoord(){return coord;}
  vector<double> getPosition(){return position;}
  vector<double> getVelocity(){return velocity;}
  vector<double> getMomentum(){return momentum;}
  vector<double> getFluidVelocity(){return fluid_velocity;}
  vector<double> getFluidThermal(){return fluid_thermal;}
  double getMass(){return masstable[id];}
  double getPathLength(){return pathlength;}
  double getQperp(){return Qperp;}

  // set functions
  void setID(int _id){id = _id;}
  void setCoord(vector<int>& _coord){
    coord = _coord;
    position = coord2position(_coord);
  }
  void setPosition(vector<double>& _position){
    position = _position;
    coord = position2coord(_position);
  }
  void setVelocity(vector<double>& _velocity);
  void setMomentum(vector<double>& _momentum);
  void setPathLength(double _pathlength){pathlength = _pathlength;}
  void addPathLength(double ds){pathlength += ds;}
  void setQperp(double _Qperp){Qperp = _Qperp;}
  void addQperp(double dQperp){Qperp += dQperp;}

  // other functions
  void stream(double dt);
  void getFluidCellData(medium* oscar, vector<cell*>* cell_array);
  void printPosition(char filename[]);
};

#endif
