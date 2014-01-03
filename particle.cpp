#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "particle.h"

particle::particle(ParameterReader* _paraRdr)
{
  paraRdr = _paraRdr;

  // particle mass table
  masstable[-4] = paraRdr->getVal("charm_mass"); // cbar
  masstable[-3] = 0.0; // sbar
  masstable[-2] = 0.0; // dbar
  masstable[-1] = 0.0; // ubar
  masstable[0] = 0.0;  // g
  masstable[1] = 0.0;  // u
  masstable[2] = 0.0;  // d 
  masstable[3] = 0.0;  // s
  masstable[4] = paraRdr->getVal("charm_mass"); // c
}

particle::~particle()
{
}

void particle::setVelocity(vector<double>& _velocity)
{
  // set velocity
  velocity = _velocity;

  // set momentum (if possible)
  vector<double> _momentum(4,0.0);
  double mass = masstable[id];
  if(mass==0.0){
    // do nothing (four-momentum cannot be determined from velocity)
  }
  else if(mass>0.0){
    for(int i=0; i<4; i++){
      _momentum[i] = (_velocity[i]*mass);
    }
  }
  else{
    cout << "error: negative mass" << endl;
  }
  momentum = _momentum;
}

void particle::setMomentum(vector<double>& _momentum)
{
  // set momentum
  momentum = _momentum;
  double pmagsq=0.0, pmag=0.0;
  for(int i=1; i<4; i++){
    pmagsq += momentum[i]*momentum[i];
  }
  pmag = pow(pmagsq,0.5);

  // set velocity (if possible)
  vector<double> _velocity(4,0.0);
  double mass = masstable[id];
  if(mass==0.0){
    _velocity[0] = 1.0; // computational trick 
    for(int i=1; i<4; i++){
      _velocity[i] = momentum[i]/pmag;
    }
  }
  else if(mass>0.0){
    for(int i=0; i<4; i++){
      _velocity[i] = (_momentum[i]/mass);
    }
  }
  else{
    cout << "error: negative mass" << endl;
  }
  velocity = _velocity;
}

void particle::stream(double dtau)
{
  double mass = masstable[id];
  vector<double> new_position(4,0.0);
  double gamma = velocity[0];
  double dt = dtau/gamma;

  for(int i=0; i<4; i++){
    new_position[i] = (position[i] + velocity[i]*clight*dt);
  }
  setPosition(new_position);
}

void particle::getFluidCellData(medium* oscar, vector<cell*>* cell_array)
{
   vector<int> coord = getCoord();
   int ifc = oscar->grabCellfromCoord(coord);
   cell* fc = (*cell_array)[ifc];
   fluid_velocity = fc->velocity;
   fluid_thermal = fc->thermal;
}

void particle::printPosition(char filename[])
{
  ofstream outfile;
  outfile.open(filename, std::ios_base::app);
  for(int i=0; i<4; i++){
    outfile << setprecision(12) << setw(22) << position[i];
  }
  outfile << endl;
  outfile.close();
}
