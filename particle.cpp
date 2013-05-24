#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include "particle.h"

particle::particle(ParameterReader* _paraRdr, pdata* _pdata)
{
  if(_pdata->coord.size() == 0) _pdata->coord = position2coord(_pdata->position);
  if(_pdata->position.size() == 0) _pdata->position = coord2position(_pdata->coord);
  paraRdr = _paraRdr;
  hquark = _pdata;
}

particle::~particle()
{
  delete hquark;
}

void particle::stream(double dtau)
{
  double gamma = hquark->velocity[0];
  double dt = dtau/gamma;

  //Update four_position
  for(int i=0; i<hquark->position.size(); i++){
    hquark->position[i] += hquark->velocity[i]*clight*dt;
  }
  //Update grid_coord
  hquark->coord = position2coord(hquark->position);
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
