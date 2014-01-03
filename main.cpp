#include <iostream>
#include <fstream> 
#include <iomanip> 
#include <stdio.h>
#include <sys/time.h>
#include "routines.h"
#include "checks.h"
using namespace std;

int main(int argc, char *argv[])
{
  // read parameters
  ParameterReader paraRdr;
  paraRdr.readFromFile("/home/morelandjs/Research/LinBolt/src/parameters.dat");
  paraRdr.readFromArguments(argc, argv);
  cout << endl; paraRdr.echo();

  // set random seed
  timeval a;
  gettimeofday(&a, 0);
  int randomSeed=a.tv_usec;
  srand48(randomSeed);

  // run routine
  int probe = 4;
  double E0=10.0, T=0.370;
  routines* energyloss  = new routines(&paraRdr);
  energyloss->transportcoeff(probe,E0,T);
  delete energyloss;
  /*checks* rate = new checks(&paraRdr);
  rate->printGamma(0);
  delete rate;*/
  

  return 1;
}

