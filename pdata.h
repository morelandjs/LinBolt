#ifndef PDATA_h
#define PDATA_h

#include <vector>
using std::vector;

struct pdata
{
  int particle_id;
  vector<int> coord;
  vector<double> position;
  vector<double> velocity;
  vector<double> momentum;
};
  
#endif
