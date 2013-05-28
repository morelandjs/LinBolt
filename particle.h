#ifndef PARTICLE_H
#define PARTICLE_H

#include "system.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "pdata.h"
#include <vector>
#include <map>
using std::vector;
using std::map;
using namespace std;

class particle : public system
{
 private:
  ParameterReader* paraRdr;
  double alphas, gs;
  map<int,double> masstable;
  
 public:
  
  pdata* pd;

  particle(ParameterReader* _paraRdr, pdata* _pdata);
  ~particle();
  void syncdata();
  vector<double> velocity2momentum(vector<double> &velocity);
  vector<double> momentum2velocity(vector<double> &momentum);
  int get_particle_id() {return pd->particle_id;}
  vector<int> get_coord() {return pd->coord;}
  vector<double> get_position() {return pd->position;}
  vector<double> get_velocity() {return pd->velocity;}
  vector<double> get_momentum() {return pd->momentum;}
  void set_particle_id(int _particle_id) {pd->particle_id=_particle_id;}
  void set_coord(vector<int> _coord) {pd->coord=_coord;}
  void set_position(vector<double> _position) {pd->position=_position;}
  void set_velocity(vector<double> _velocity) {pd->velocity=_velocity;}
  void set_momentum(vector<double> _momentum) {pd->momentum=_momentum;}
  void stream(double dt);
  void printPosition(char filename[]);
};

#endif
