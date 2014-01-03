#include <iostream>
#include <fstream> 
#include <iomanip> 
#include <stdio.h>
#include <sys/time.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include "routines.h"
using std::vector;
using namespace std;


routines::routines(ParameterReader* _paraRdr)
{
  // fix parameters
  paraRdr = _paraRdr;
  verbose = paraRdr->getVal("verbose");
  mc = paraRdr->getVal("charm_mass");
  alphas = paraRdr->getVal("alphas");
  gs = pow(4.0*M_PI*alphas,0.5);
  tmax = paraRdr->getVal("grid_t1");
  tmin = paraRdr->getVal("grid_t0");
  nt = paraRdr->getVal("grid_nt");
  dtau = (tmax-tmin)/double(nt);
}

void routines::transportcoeff(int species, double E0, double T)
{
  int nparticles = 1000, bins = 30;
  double xmin = 0, xmax = 30.0, dx = (xmax-xmin)/double(bins); 
  double energy = 0.0, pl = 0.0, ptsq = 0.0, dptsq = 0.0, pathlength = 0.0, length = 0.0, time = 0.0;
  vector<double> vi(4,0.0), vf(4,0.0), pi(4,0.0), pf(4,0.0);
  double vin, dt, pin, pfn;

  // initialize classes
  medium* oscar = new medium(paraRdr, "/run/media/morelandjs/b09c7850-22ca-4fc5-b7a3-96981704d524/morelandjs/OSCAR2008H.dat");
  scattering* dynamics = new scattering(paraRdr);

  // create probes
  for(int ispecies=0; ispecies<nparticles; ispecies++){
    particle* probe = new particle(paraRdr); probe->setID(species);
    vector<double> position(4,0.0), momentum(4,0.0);
    double mass = probe->getMass();
    double Emag = E0, Pmag = sqrt(Emag*Emag-mass*mass);
    position[0] = 0.0; position[1] = 0.0; position[2] = 0.0; position[3] = 0.0;
    momentum[0] = Emag; momentum[1] = Pmag; momentum[2] = 0.0; momentum[3] = 0.0;
    probe->setPosition(position);
    probe->setMomentum(momentum);
    probe->setPathLength(0.0);
    probe->setQperp(0.0);
    particle_list.push_back(probe);
  }

  // open output file
  ofstream energyloss, qhat, pcoord;
  energyloss.open("/home/morelandjs/Research/LinBolt/results/energyloss.dat", std::ios_base::trunc);
  qhat.open("/home/morelandjs/Research/LinBolt/results/qhat.dat", std::ios_base::trunc);
  pcoord.open("/home/morelandjs/Research/LinBolt/results/particle.dat", std::ios_base::trunc);

  // populate the static medium
  oscar->populateStatic();
  cell_array = oscar->getHydro();

  // create histogram bins
  vector<double> energy_bins(bins,0.0);
  vector<double> ptsq_bins(bins,0.0);
  vector<int> ebin_counts(bins,0);
  vector<int> pbin_counts(bins,0);
  vector<double> distance_bins(bins,0.0);
  for(int i=0; i<bins; i++){
    distance_bins[i] = xmin+i*dx+dx/2.0;
  }

  cout << "Starting time evolution:" << endl;

  // time loop
  for(int itime=0; itime<nt; itime++){
    //particle loop
    for(int n=0; n<particle_list.size(); n++){

      // free stream 
      pn = particle_list[n];
      vi = pn->getVelocity();
      pi = pn->getMomentum();
      pin = pow(pow(pi[1],2.0)+pow(pi[2],2.0)+pow(pi[3],2.0),0.5);
      dt = dtau/vi[0];
      vin = pow(pow(vi[1],2.0)+pow(vi[2],2.0)+pow(vi[3],2.0),0.5);
      pn->addPathLength(vin*dt);
      pn->stream(dtau);

      // get energy and distance traveled
      energy = pn->getMomentum()[0];
      ptsq = pn->getQperp();
      pathlength = pn->getPathLength();
      time = pn->getPosition()[0];
      length = pn->getPosition()[1];

      // bin particles for energy loss and qhat
      for(int i=0; i<bins; i++){
	if(distance_bins[i] <= pathlength && pathlength < distance_bins[i+1]){
	  energy_bins[i] += energy;
	  ebin_counts[i] ++;
	}
	if(distance_bins[i] <= pathlength && pathlength < distance_bins[i+1]){
	  ptsq_bins[i] += ptsq/pathlength;
	  pbin_counts[i] ++;
	}
      }

      if(itime%5 == 0 && pn->getPosition()[1]<50.0 && abs(pn->getPosition()[2])<15.0){
      pcoord << setprecision(12) << setw(22) << pn->getPosition()[0] << setprecision(12) << setw(22) << pn->getPosition()[1] << setprecision(12) 
             << setw(22) << pn->getPosition()[2] << setprecision(12) << setw(22) << pn->getPosition()[3] << endl;
      }

      // scatter
      pn->getFluidCellData(oscar,cell_array);
      if(dynamics->isWounded(pn,dtau)){
	vector<particle*> products = dynamics->scatter(pn);
	int id3 = products[0]->getID();
	int id4 = products[1]->getID();
	double E3 = products[0]->getMomentum()[0];
	double E4 = products[1]->getMomentum()[0];
	if((id3==species && id4==species && E3>=E4) || (id3==species && id4 != species)){
	  pf = products[0]->getMomentum();
	  pfn = pow(pow(pf[1],2.0)+pow(pf[2],2.0)+pow(pf[3],2.0),0.5);
	  pl = (pi[1]*pf[1]+pi[2]*pf[2]+pi[3]*pf[3])/pin;
	  dptsq = pow(pfn,2.0)-pow(pl,2.0);
	  products[0]->setQperp(ptsq);
	  products[0]->addQperp(dptsq);
	  products[0]->setPathLength(pathlength);
	  particle_list_buffer.push_back(products[0]);
	}
	else if((id3==species && id4==species && E3<E4) || (id4==species && id3 != species)){
	  pf = products[1]->getMomentum();
	  pfn = pow(pow(pf[1],2.0)+pow(pf[2],2.0)+pow(pf[3],2.0),0.5);
	  pl = (pi[1]*pf[1]+pi[2]*pf[2]+pi[3]*pf[3])/pin;
	  dptsq = pow(pfn,2.0)-pow(pl,2.0);
	  products[1]->setQperp(ptsq);
	  products[1]->addQperp(dptsq);
	  products[1]->setPathLength(pathlength);
	  particle_list_buffer.push_back(products[1]);
	}
      }
      else particle_list_buffer.push_back(pn);

    }

    // transfer new particles to the buffer
    particle_list = particle_list_buffer; 
    particle_list_buffer.clear();
    drawProgressBar(60,double(itime)/double(nt));
  }

  for(int i=0; i<bins; i++){
      energyloss << setprecision(12) << setw(22) << distance_bins[i] << setprecision(12) << setw(22) << energy_bins[i]/double(ebin_counts[i]) << endl;
      qhat << setprecision(12) << setw(22) << distance_bins[i] << setprecision(12) << setw(22) << ptsq_bins[i]/double(pbin_counts[i]) << endl;
    }

  // clean up
  energyloss.close();
  qhat.close();
  pcoord.close();
  particle_list.clear();
  particle_list_buffer.clear();
  delete dynamics;
  delete oscar;
}

void routines::charminhydro(double E0)
{
  int nparticles = 100;

  // Initialize classes
  medium* oscar = new medium(paraRdr, "/run/media/morelandjs/b09c7850-22ca-4fc5-b7a3-96981704d524/morelandjs/OSCAR2008H.dat");
  scattering* dynamics = new scattering(paraRdr);

  for(int ispecies=0; ispecies<nparticles; ispecies++){
    particle* probe = new particle(paraRdr); probe->setID(4);
    vector<double> position(4,0.0), momentum(4,0.0);
    double Emag = E0, Pmag = sqrt(Emag*Emag-mc*mc);
    position[0] = 0.0; position[1] = -30.0; position[2] = 0.0; position[3] = 0.0;
    momentum[0] = Emag; momentum[1] = Pmag; momentum[2] = 0.0; momentum[3] = 0.0;
    probe->setPosition(position);
    probe->setMomentum(momentum);
    particle_list.push_back(probe);
  }

  //dynamics->print_rates(); // test 2->2 rates
  
  // OSCAR time loop
  //========================================================
  for(int imin=0; imin<300; imin+=20){

    oscar->populateOSCAR(imin,imin+20);
    cell_array = oscar->getHydro();
    
    for(int itime=imin; itime<(imin+20); itime+=1){
      for(int n=0; n<particle_list.size(); n++){

	pn = particle_list[n];
	pn->stream(dtau);
	pn->getFluidCellData(oscar,cell_array);

	if(dynamics->isWounded(pn,dtau)){
	  if(verbose) cout << "** " << particle_list.size() << " particles **" << endl;
	  vector<particle*> products = dynamics->scatter(pn);
	  if(products[0]->getID()==4) particle_list_buffer.push_back(products[0]);
	  else if(products[1]->getID()==4) particle_list_buffer.push_back(products[1]);
	}
	else particle_list_buffer.push_back(pn);
      }

      particle_list = particle_list_buffer;
      particle_list_buffer.clear();

      oscar->populate2D(itime,0);
      sprintf(path,"/home/morelandjs/morelandjs/LinBolt/data/energy_dens_%05d.dat",itime); 
      oscar->print2D(path);
    }
  }
}
