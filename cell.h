#ifndef CELL_h
#define CELL_h

#include <vector>
using std::vector;

struct cell
{
  vector<int> coord;
  vector<double> position;
  vector<double> lagrange;
  vector<double> velocity;
  vector<double> thermal;
};
  
#endif
