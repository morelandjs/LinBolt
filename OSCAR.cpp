#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include "OSCAR.h"
#include "arsenal.h"

OSCAR::OSCAR(ParameterReader *paraRdr_in, char _filename[])
{
  //Initialize parameter reader
  paraRdr = paraRdr_in;

  //Initialize OSCAR file stream
  string line;
  filename = _filename;
  ifstream infile;   
  infile.open(filename);

  //Skip the OSCAR file header data
  if (!infile)
    cout << "Error opening data file." << endl;  
  else{
    for(int l=0; l<14; l++){getline(infile,line);}
  }
  infile.close();

  // Read coord system parameters from OSCAR header
  vector<int> _dim;
  vector<double> _gridmin, _gridmax;

  nt = paraRdr->getVal("grid_nt");
  _dim.push_back(nt);
  t0 = paraRdr->getVal("grid_t0");
  _gridmin.push_back(t0);
  t1 = paraRdr->getVal("grid_t1");
  _gridmax.push_back(t1);
  nx = paraRdr->getVal("grid_nx");
  _dim.push_back(nx);
  x0 = paraRdr->getVal("grid_x0");
  _gridmin.push_back(x0);
  x1 = paraRdr->getVal("grid_x1");
  _gridmax.push_back(x1);
  ny = paraRdr->getVal("grid_ny");
  _dim.push_back(ny);
  y0 = paraRdr->getVal("grid_y0");
  _gridmin.push_back(y0);
  y1 = paraRdr->getVal("grid_y1");
  _gridmax.push_back(y1);
  nz = paraRdr->getVal("grid_nz");
  _dim.push_back(nz);
  z0 = paraRdr->getVal("grid_z0");
  _gridmin.push_back(z0);
  z1 = paraRdr->getVal("grid_z1");
  _gridmax.push_back(z1);
  CC = paraRdr->getVal("dim_CC");
  DD = paraRdr->getVal("dim_DD");
  TT = paraRdr->getVal("dim_TT");

  // Fix the system parameters
  setDim(_dim);
  setGridMin(_gridmin);
  setGridMax(_gridmax);

  // Read OSCAR parameter file
  parse_it = paraRdr->getVal("parse_it");
  parse_ix = paraRdr->getVal("parse_ix");
  parse_iy = paraRdr->getVal("parse_iy");
  parse_iz = paraRdr->getVal("parse_iz");
  parse_tau = paraRdr->getVal("parse_tau");
  parse_x = paraRdr->getVal("parse_x");
  parse_y = paraRdr->getVal("parse_y");
  parse_eta = paraRdr->getVal("parse_eta");
  parse_e = paraRdr->getVal("parse_e");
  parse_p = paraRdr->getVal("parse_p");
  parse_T = paraRdr->getVal("parse_T");
  parse_R_qgp = paraRdr->getVal("parse_R_qgp");
  parse_vx = paraRdr->getVal("parse_vx");
  parse_vy = paraRdr->getVal("parse_vy");
  parse_y_L = paraRdr->getVal("parse_y_L");
  parse_n = paraRdr->getVal("parse_n");
  parse_mu = paraRdr->getVal("parse_mu");
  parse_Diss = paraRdr->getVal("parse_Diss");
  parse_Tr = paraRdr->getVal("parse_Tr");
  toggle_it = paraRdr->getVal("toggle_it");
  toggle_ix = paraRdr->getVal("toggle_ix");
  toggle_iy = paraRdr->getVal("toggle_iy");
  toggle_iz = paraRdr->getVal("toggle_iz");
  toggle_tau = paraRdr->getVal("toggle_tau");
  toggle_x = paraRdr->getVal("toggle_x");
  toggle_y = paraRdr->getVal("toggle_y");
  toggle_eta = paraRdr->getVal("toggle_eta");
  toggle_e = paraRdr->getVal("toggle_e");
  toggle_p = paraRdr->getVal("toggle_p");
  toggle_T = paraRdr->getVal("toggle_T");
  toggle_R_qgp = paraRdr->getVal("toggle_R_qgp");
  toggle_vx = paraRdr->getVal("toggle_vx");
  toggle_vy = paraRdr->getVal("toggle_vy");
  toggle_y_L = paraRdr->getVal("toggle_y_L");
  toggle_n = paraRdr->getVal("toggle_n");
  toggle_mu = paraRdr->getVal("toggle_mu");
  toggle_Diss = paraRdr->getVal("toggle_Diss");
  toggle_Tr = paraRdr->getVal("toggle_Tr");

  // Initialize the 2D array saver
  array2D = new double*[nx];
  for(int ix=0;ix<nx;ix++){
    array2D[ix] = new double[ny]();
    for (int iy=0;iy<ny;iy++) array2D[ix][iy]=0.0;
  }
}

OSCAR::~OSCAR()
{
  for(int iy=0;iy<ny;iy++) delete [] array2D[iy];
  delete [] array2D;

  //Deconstruct OSCAR array
  hydro.clear();
}

void OSCAR::drawProgressBar(int len, double percent)
{
  cout << "\x1B[2K"; // Erase the entire current line.
  cout << "\x1B[0E"; // Move to the beginning of the current line.
  string progress;
  for (int i = 0; i < len; ++i) {
    if (i < static_cast<int>(len*percent+0.5)) {
      progress += "=";
    } else {
      progress += " ";
    }
  }
  cout << "[" << progress << "] " << (static_cast<int>(100*percent+0.5)) << "%";
  flush(cout); // Required.
}

void OSCAR::populateOSCAR(ParameterReader *paraRdr_in, int _ntmin, int _ntmax)
{
  cout << endl;
  cout << "populating the hydro array: " << "itime " << _ntmin << " to " << _ntmax << endl;
  hydro.clear();
  ntmin = _ntmin;
  ntmax = _ntmax;

  //Initialize parameter reader
  paraRdr = paraRdr_in;
  
  //Open the file stream
  string line;
  ifstream infile;   
  infile.open(filename);

  //Skip the header data
  if (!infile)
    cout << "Error opening data file." << endl;  
  else{
    for(int l=0; l<14; l++) getline(infile,line);
  }

  //Skip to the initial read position
  int start_read = ntmin*nx*ny*nz;
  int reads = (ntmax-ntmin+1)*nx*ny*nz;
  for(int i=0; i<start_read;i++) getline(infile,line);

  //Read the OSCAR file
  for(int i=0; i<reads; i++)
    {
      double perc = (double)i/(double)reads;
      int q = reads/40;
      if(i%q==0)drawProgressBar(40,perc);
      // Read in a hydro cell
      if(parse_it) {infile >> it;}
      if(parse_ix) {infile >> ix;}
      if(parse_iy) {infile >> iy;}
      if(parse_iz) {infile >> iz;}
      if(parse_tau) {infile >> tau;}
      if(parse_x) {infile >> x;}
      if(parse_y) {infile >> y;}
      if(parse_eta) {infile >> eta;}
      if(parse_e) {infile >> e;}
      if(parse_p) {infile >> p;}
      if(parse_T) {infile >> T;}
      if(parse_R_qgp) {infile >> R_qgp;}
      if(parse_vx) {infile >> vx;}
      if(parse_vy) {infile >> vy;}
      if(parse_y_L) {infile >> y_L;}
      infile.ignore(1000, '\n');

      cell* fluidcell = new cell;
      // Store the spacetime indices and coords
      if(toggle_it) {fluidcell->coord.push_back(it);}
      if(toggle_ix) {fluidcell->coord.push_back(ix);}
      if(toggle_iy) {fluidcell->coord.push_back(iy);}
      if(toggle_iz) {fluidcell->coord.push_back(iz);}
      fluidcell->position = coord2position(fluidcell->coord);
	  
      // Store additional Lagrangian coordinates
      if(toggle_tau) fluidcell->lagrange.push_back(tau);
      if(toggle_x) fluidcell->lagrange.push_back(x);
      if(toggle_y) fluidcell->lagrange.push_back(y);
      if(toggle_eta) fluidcell->lagrange.push_back(eta);

      // Store thermal properties
      if(toggle_e) fluidcell->thermal.push_back(e);
      if(toggle_p) fluidcell->thermal.push_back(p);
      if(toggle_T) fluidcell->thermal.push_back(T);
      if(toggle_R_qgp) fluidcell->thermal.push_back(R_qgp);

      // Store flow velocities
      double gamma_perp=1;
      if(toggle_vx && toggle_vy && toggle_y_L) gamma_perp= 1/pow(1-vx*vx-vy*vy,0.5);
      if(toggle_vx && toggle_vy && toggle_y_L) fluidcell->velocity.push_back(gamma_perp*cosh(y_L));
      if(toggle_vx && toggle_vy && toggle_y_L) fluidcell->velocity.push_back(gamma_perp*vx);
      if(toggle_vx && toggle_vy && toggle_y_L) fluidcell->velocity.push_back(gamma_perp*vy);
      if(toggle_vx && toggle_vy && toggle_y_L) fluidcell->velocity.push_back(gamma_perp*sinh(y_L));

      // Add the fluid cell to the hydro vector
      hydro.push_back(fluidcell);
    }
  cout << endl;
  infile.close();
}

// Find the vector index for the given hydro cell
int OSCAR::grabCellfromCoord(vector<int> _coord)
{
  int ncell=0;
  if(toggle_it==true) ncell += (_coord[0]-ntmin)*nx*ny*nz;
  if(toggle_ix==true) ncell += _coord[1]*ny*nz;
  if(toggle_iy==true) ncell += _coord[2]*nz;
  if(toggle_iz==true) ncell += _coord[3];
  return ncell; 
}

void OSCAR::populate2D(int _it, int _thermalID)
{
  int ncell;
  vector<int> coord;
  
  for(int ix=0; ix<nx; ix++){
    for(int iy=0; iy<ny; iy++){
      coord.clear();
      coord.push_back(_it);
      coord.push_back(ix);
      coord.push_back(iy);
      ncell = this->grabCellfromCoord(coord);
      array2D[ix][iy] = hydro[ncell]->thermal[_thermalID];
    }
  }
}

void OSCAR::print2D(char filename[])
{
  vector<int> coord;
  ofstream outfile;
  outfile.open(filename, std::ios_base::trunc);

  for(int ix=0;ix<nx;ix++){
    for(int iy=0;iy<ny;iy++){
      coord.clear();
      coord.push_back(0);
      coord.push_back(ix);
      coord.push_back(iy);
      vector<double> position = coord2position(coord);
      outfile << setprecision(12) << setw(22) << position[1] 
	      << setprecision(12) << setw(22) << position[2]
	      << setprecision(12) << setw(22) << array2D[ix][iy] << endl;
    }
    outfile << endl;
  }
  outfile.close();
}

vector<double> OSCAR::sample_boltzmann(double mass, double temp)
{
  vector<double> p(4, 0.0);
  double a,b,c;
  double r0,r1,r2,r3;
  double I1,I2,I3,Itot;
  double pmag,E,ctheta,stheta,phi,pgauss,K,f;

  if(temp/mass>0.6){
    cout << 1 << endl;
    do{
      cout << 2 << endl;
      r1=drand48();
      r2=drand48();
      r3=drand48();
      a=-log(r1); b=-log(r2); c=-log(r3);
      cout << "a: " << a << endl;
      cout << "b: " << b << endl;
      cout << "c: " << c << endl;
      pmag=temp*(a+b+c);
      cout << "pmag: " << pmag << endl;
      cout << "mass: " << mass << endl;
      cout << pow(pmag*pmag+mass*mass,0.5) << endl;
      p[0]=pow(pmag*pmag+mass*mass,0.5);
      cout << "p0: " << p[0] << endl;
      cout << "temp: " << temp << endl;
      f = 1/(exp((p[0]-pmag)/temp)+1);
      cout << "f: " << f << endl;
    } while(drand48()>f);
    cout << 3 << endl;
    ctheta=(a-b)/(a+b);
    stheta=sqrt(1.0-ctheta*ctheta);
    phi=temp*temp*pow(a+b,2)/(pmag*pmag);
    phi=2.0*M_PI*phi;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
    p[3]=pmag*ctheta; 
  }
  else{
    cout << 4 << endl;
    do{
      cout << 5 << endl;
      r0=drand48();
      I1=mass*mass;
      I2=2.0*mass*temp;
      I3=2.0*temp*temp;
      Itot=I1+I2+I3;
      if(r0<I1/Itot){
	cout << 6 << endl;
	r1=drand48();
	K=-temp*log(r1);
      }
      else if(r0<(I1+I2)/Itot){
	cout << 7 << endl;
	r1=drand48();
	r2=drand48();
	K=-temp*log(r1*r2);
      }
      else{
	cout << 8 << endl;
	r1=drand48();
	r2=drand48();
	r3=drand48();
	K=-temp*log(r1*r2*r3);
      }
      cout << 9 << endl;
      E=K+mass;
      pmag=sqrt(E*E-mass*mass);
      r0=drand48();
    } while(r0>pmag/E);
    cout << 10 << endl;
    phi=2.0*M_PI*drand48();
    cout << "phi: " << phi << endl;
    ctheta=1.0-2.0*drand48();
    cout << "ctheta: " << ctheta << endl;
    stheta=sqrt(1.0-ctheta*ctheta);
    cout << "stheta: " << stheta << endl;
    cout << "E: " << E << endl;
    cout << "pmag: " << pmag << endl;
    p[0]=E;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
    p[3]=pmag*ctheta;

    cout << E << endl;
    cout << pmag << " " << stheta << " " << cos(phi) << endl;
    cout << pmag*stheta*cos(phi) << endl;
    cout << pmag*stheta*sin(phi) << endl;
    cout << pmag*ctheta << endl;
  }
  return p;
}








