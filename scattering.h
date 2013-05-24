#ifndef SCATTERING_H
#define SCATTERING_H

#include "system.h"
#include "ParameterReader.h"
#include "particle.h"
#include "OSCAR.h"
using namespace std;

class scattering : public system
{
 private:
  ParameterReader* paraRdr;
 public:
  scattering(ParameterReader* _paraRdr);
  ~scattering();
  double sampleDiffXS(double s, double temp, int i, int j, int k, int l);
  double invertDiffXS(double s, double temp, int i, int j, int k, int l, double IM2);
  double Gamma_i(double E1, double temp, int i);
  double Gamma_ij2kl(double E1, double temp, int i, int j, int k, int l);
  double Sigma_ij2kl(double E1, double E2, double m, int i, int j, int k, int l);
  double M2_ij2kl_integrated(double s, double t, int i, int j, int k, int l);
  double M2_ij2kl(double s, double t, double u, int i, int j, int k, int l);
  double F(double E2, double temp);
  double A10(double E1, double E2, double m);
  double A11(double E1, double E2, double m);
  double A12(double E1, double E2, double m);
  double A13(double E1, double E2, double m);
  double A14(double E1, double E2, double m);
  double A15(double E1, double E2, double m);
  double A16(double E1, double E2, double m);
  double A17(double E1, double E2, double m);
  double A18(double E1, double E2, double m);
  void printGamma(char filename[]);
 protected:
  double alphas, gs, sm, gq;
  vector<double> CDF;
};

#endif
