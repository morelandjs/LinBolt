#ifndef SCATTERING_H
#define SCATTERING_H

#include "system.h"
#include "arsenal.h"
#include "medium.h"
#include "particle.h"
#include "ParameterReader.h"
#include <vector>
using std::vector;
using namespace std;

class scattering : public system
{
 private:
  ParameterReader* paraRdr;
  bool verbose, gluon, quark, hquark;
  double E1_min, E1_max, dE1, E2_min, E2_max, dE2, T_min, T_max, dT;
  int nE1, nE2, nT, ns_table;

 public:
  double **Gamma_gggg, **Gamma_ggqqbar, **Gamma_gqgq, **Gamma_qiqjqiqj,**Gamma_qiqiqiqi, 
    **Gamma_qiqibarqjqjbar, **Gamma_qiqibarqiqibar, **Gamma_qqbargg, ** Gamma_qgqg, 
    **Gamma_cqcq, **Gamma_cgcg;

  double **Omega_gggg, **Omega_ggqqbar, **Omega_gqgq, **Omega_qiqjqiqj, **Omega_qiqiqiqi, 
    **Omega_qiqibarqjqjbar, **Omega_qiqibarqiqibar, **Omega_qqbargg, **Omega_qgqg, 
    **Omega_cqcq, **Omega_cgcg;

  scattering(ParameterReader* _paraRdr);
  ~scattering();

  // Collision Determination Routines
  bool isWounded(particle* p1, double dtau);
  vector<particle*> scatter(particle* p1);

  // Sampling Routines (energy, angles, etc)
  vector<int> sample_2to2(double E1, double T, int i);
  double sample_E2(double E1, double temp, vector<int> ij2kl);
  double dGammadE2(double E1, double E2, double temp, vector<int> ij2kl);
  double sample_theta12(double E1, double E2, double temp, vector<int> ij2kl);
  double sample_theta13(double s, double pmag, double T, vector<int> ij2kl);
  double sample_phi12();
  double sample_phi13();
  
  // Gamma Function Routines
  double incl_Gamma(double E1, double T, int i);
  double Gamma(double E1, double T, vector<int> ij2kl);
  double interp_Gamma(double E1, double T, double **Gamma);
  double** populate_Gamma(vector<int> ij2kl);
  double calc_Gamma(double E1, double T, vector<int> ij2kl);
 
  // sigma Function Routines
  double invert_sigmaI(double s, double T, vector<int> ij2kl, double intMsq);
  double sigmaI(double s, double t, double temp, vector<int> ij2kl);

  // Omega Function Routines
  double** populate_OmegaI(vector<int> ij2kl);
  double invert_OmegaI(double E1, double E2, double temp, vector<int> ij2kl, double frac);
  double OmegaI(double s, double temp, vector<int> ij2kl);
  double Omega(double E1, double E2, double temp, vector<int> ij2kl);
  double interp_OmegaI(double s, double temp, vector<int> ij2kl, double **OmegaI);
  double calc_OmegaI(double s, double temp, vector<int> ij2kl);

  // Kinematic Limits
  double get_cutoff(double T);
  vector<double> get_s_limits(double E1, double E2, double Q0, vector<int> ij2kl);
  vector<double> get_t_limits(double s, double Q0, vector<int> ij2kl);

  // Mandelstam & Dist. Functions
  double F(double E2, double T, int j);
  double get_MandelstamS(vector<double> &p1, vector<double> &p2);
  double get_MandelstamT(vector<double> &p1, vector<double> &p3);
  double get_MandelstamU(vector<double> &p1, vector<double> &p4);

  // Global Constants
 protected:
  double alphas, gs, gg, gq, mc, db;
};

#endif
