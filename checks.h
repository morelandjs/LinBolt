#ifndef CHECKS_H
#define CHECKS_H

#include "system.h"
#include "medium.h"
#include "particle.h"
#include "scattering.h"
#include "arsenal.h"
#include "ParameterReader.h"

class checks
{
 private:
  int nE1, nT;
  double E1_min, E1_max, dE1, T_min, T_max, dT;
  ParameterReader* paraRdr;
  scattering* dynamics;
  
 public:
  checks(ParameterReader* _paraRdr);
  void printGamma(int species);
};

#endif
