#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "particle.h"

particle::particle(ParameterReader* _paraRdr, pdata* _pdata)
{
  paraRdr = _paraRdr;
  pd = _pdata;
  // particle mass table
  masstable[-2] = 0.0048; // dbar
  masstable[-1] = 0.0024; // ubar
  masstable[0] = 0.0;     // g
  masstable[1] = 0.0024;  // u
  masstable[2] = 0.0048;  // d 
  syncdata();
}

particle::~particle()
{
}

void particle::syncdata()
{
  int cs = pd->coord.size();
  int xs = pd->position.size();
  int vs = pd->velocity.size();
  int ps = pd->momentum.size();
  if(cs == 0 && xs != 0) set_coord(position2coord(pd->position));
  if(xs == 0 && cs != 0) set_position(coord2position(pd->coord));
  if(vs == 0 && ps != 0) set_velocity(momentum2velocity(pd->momentum));
  if(ps == 0 && vs != 0) set_momentum(velocity2momentum(pd->velocity));
}

vector<double> particle::velocity2momentum(vector<double> &velocity)
{
  double mass = masstable[pd->particle_id];
  vector<double> momentum;
  int vs = velocity.size();
  
  if(mass==0.0){
    // do nothing (massless four-velocity is undefined)
  }
  else if(mass>0.0){
    for(int i=0; i<vs; i++){
      momentum.push_back(velocity[i]*mass);
    }
  }
  else{
    cout << "error: negative mass" << endl;
  }
  return momentum; 
}

vector<double> particle::momentum2velocity(vector<double> &momentum)
{
  double mass = masstable[pd->particle_id];
  vector<double> velocity;
  int ps = momentum.size();

  if(mass==0.0){
    // do nothing (massless four-velocity is undefined)
  }
  else if(mass>0.0){
    for(int i=0; i<ps; i++){
      velocity.push_back(momentum[i]/mass);
    }
  }
  else{
    cout << "error: negative mass" << endl;
  }
  return velocity; 
}

void particle::stream(double dtau)
{
  syncdata();
  double mass = masstable[pd->particle_id];

  if(mass==0.0){
    for(int i=0; i<pd->position.size(); i++){
      pd->position[i] += dtau;
    }
  }
  else if(mass>0.0){
    double gamma = pd->velocity[0];
    double dt = dtau/gamma;
    for(int i=0; i<pd->position.size(); i++){
      pd->position[i] += pd->velocity[i]*clight*dt;
    }
  }
  else{
    cout << "error: negative mass" << endl;
  }

  syncdata();
}

void particle::printPosition(char filename[])
{
  ofstream outfile;
  outfile.open(filename, std::ios_base::app);

  double xx = get_position()[1];
  double yy = get_position()[2];
  outfile << setprecision(12) << setw(22) << xx
	  << setprecision(12) << setw(22) << yy << endl;
  outfile.close();
}
