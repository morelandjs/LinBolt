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
  vector<int> sample2to2Processes(double s, double temp, int i);
  double Gamma_i(double E1, double temp, int i);
  double Gamma_ij2kl(double E1, double temp, int i, int j, int k, int l);
  double Sigma_ij2kl(double E1, double E2, double m, int i, int j, int k, int l);
  double samplePhi();
  double sampleTheta(double s, double temp, int i, int j, int k, int l);
  double invert_IM2_ij2kl(double s, double temp, int i, int j, int k, int l, double IM2);
  double IM2_ij2kl(double s, double t, int i, int j, int k, int l);
  double M2_ij2kl(double s, double t, double u, int i, int j, int k, int l);
  double A10(double E1, double E2, double m);
  double A11(double E1, double E2, double m);
  double A12(double E1, double E2, double m);
  double A13(double E1, double E2, double m);
  double A14(double E1, double E2, double m);
  double A15(double E1, double E2, double m);
  double A16(double E1, double E2, double m);
  double A17(double E1, double E2, double m);
  double A18(double E1, double E2, double m);
  double F(double E2, double temp);
  void printGamma(char filename[]);
 protected:
  double alphas, gs, sm, gq;
};

#endif
