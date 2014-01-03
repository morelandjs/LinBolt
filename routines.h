#ifndef ROUTINES_H
#define ROUTINES_H

#include "system.h"
#include "medium.h"
#include "particle.h"
#include "scattering.h"
#include "arsenal.h"
#include "ParameterReader.h"

class routines
{
private:
  bool verbose;
  ParameterReader* paraRdr;
  double alphas, gs, mc, tmax, tmin, nt, dtau;
  
  scattering* dynamics;
  vector<cell*>* cell_array;
  vector<particle*> particle_list, particle_list_buffer;
  particle* pn;
  char path[100];
  
public:
  routines(ParameterReader* _paraRdr);
  void transportcoeff(int species, double E0, double T);
  void charminhydro(double E0);
};

#endif
