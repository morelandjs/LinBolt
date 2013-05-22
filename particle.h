#ifndef PARTICLE_H
#define PARTICLE_H

#include "system.h"
#include "ParameterReader.h"
#include "pdata.h"
#include <vector>
using std::vector;
using namespace std;

class particle : public system
{
 private:
  ParameterReader* paraRdr;
  double alphas, gs;
 public:
  pdata* hquark;
  particle(ParameterReader* _paraRdr, pdata* _pdata);
  ~particle();
  vector<int> get_particle_id() {return hquark->particle_id;}
  vector<int> get_coord() {return hquark->coord;}
  vector<double> get_position() {return hquark->position;}
  vector<double> get_velocity() {return hquark->velocity;}
  vector<double> get_momentum() {return hquark->momentum;}
  void set_particle_id(vector<int> _particle_id) {hquark->particle_id=_particle_id;}
  void set_coord(vector<int> _coord) {hquark->coord=_coord;}
  void set_position(vector<double> _position) {hquark->position=_position;}
  void set_velocity(vector<double> _velocity) {hquark->velocity=_velocity;}
  void set_momentum(vector<double> _momentum) {hquark->momentum=_momentum;}
  void stream(double dt);
  void printPosition(char filename[]);
};

#endif
