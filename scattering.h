#ifndef SCATTERING_H
#define SCATTERING_H

#include "system.h"
#include "OSCAR.h"
#include "particle.h"
#include "ParameterReader.h"
#include <vector>
using std::vector;
using namespace std;

class scattering : public system
{
 private:
  ParameterReader* paraRdr;
 public:
  scattering(ParameterReader* _paraRdr);
  ~scattering();
  vector<int> sample2to2(double s, double temp, int i);
  double Gamma_i(double E1, double temp, int i);
  double Gamma_ij2kl(double E1, double temp, vector<int> ij2kl);
  double Sigma_ij2kl(double E1, double E2, double m, vector<int> ij2kl);
  double samplePhi();
  double sampleTheta(double s, double temp, vector<int> ij2kl);
  double invert_IM2_ij2kl(double s, double temp, vector<int> ij2kl, double IM2);
  double IM2_ij2kl(double s, double t, vector<int> ij2kl);
  double M2_ij2kl(double s, double t, double u, vector<int> ij2kl);
  double A10(double E1, double E2, double m);
  double A11(double E1, double E2, double m);
  double A12(double E1, double E2, double m);
  double A13(double E1, double E2, double m);
  double A14(double E1, double E2, double m);
  double A15(double E1, double E2, double m);
  double A16(double E1, double E2, double m);
  double A17(double E1, double E2, double m);
  double A18(double E1, double E2, double m);
  double getMandelstamS(vector<double> &p1, vector<double> &p2);
  double F(double E2, double temp);
  void printGamma(char filename[]);
 protected:
  double alphas, gs, sm, gq;
};

#endif
