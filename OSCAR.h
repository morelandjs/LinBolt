#ifndef OSCAR_H
#define OSCAR_H

#include "system.h"
#include "arsenal.h"
#include "ParameterReader.h"
#include "cell.h"
#include <vector>
using std::vector;
using namespace std;

class OSCAR : public system
{
 private:
  char* filename;
  int nt, nx, ny, nz, CC, DD, TT;
  int it, ix, iy, iz;
  int ntmin, ntmax;
  double t0, t1, x0, x1, y0, y1, z0, z1;
  double tau, x, y, eta;
  double vx, vy, y_L;
  double e, p, T, R_qgp;
  bool parse_it, parse_ix, parse_iy, parse_iz, parse_tau, parse_x, parse_y;
  bool parse_eta, parse_e, parse_p, parse_T, parse_R_qgp, parse_vx, parse_vy;
  bool parse_y_L, parse_n, parse_mu, parse_Diss, parse_Tr;
  bool toggle_it, toggle_ix, toggle_iy, toggle_iz, toggle_tau, toggle_x, toggle_y;
  bool toggle_eta, toggle_e, toggle_p, toggle_T, toggle_R_qgp, toggle_vx, toggle_vy;
  bool toggle_y_L, toggle_n, toggle_mu, toggle_Diss, toggle_Tr;
  double** array2D;
  ParameterReader* paraRdr;
  vector<cell*> hydro;

 public:
  OSCAR(ParameterReader* paraRdr_in, char filename[]);
  ~OSCAR();
  void populateOSCAR(ParameterReader* paraRdr_in, int _ntmin, int _ntmax); 
  int grabCellfromCoord(vector<int> _coord);
  void populate2D(int _it, int _thermalID);
  void print2D(char filename[]);
  vector<cell*>* getHydro() {return &hydro;}
  void drawProgressBar(int len, double percent);
  vector<double> sample_boltzmann(double mass, double T);
};

#endif
