#include <iostream>
#include <stdio.h>
#include "system.h"
#include "OSCAR.h"
#include "particle.h"
#include "scattering.h"
#include "arsenal.h"
#include <cmath>
#include <vector>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
using std::vector;
using namespace boost::assign;
using namespace std;


int main(int argc, char *argv[])
{
  // vector<int> loc;
  int ithermal=0;
  char path[100];

  // Read-in parameters
  ParameterReader paraRdr;
  paraRdr.readFromFile("parameters.dat");
  paraRdr.readFromArguments(argc, argv);
  paraRdr.echo();

  // Construct the hydro object
  OSCAR* oscar = new OSCAR(&paraRdr, "/run/media/jsm55/Seagate Backup Plus Drive/Data/OSCAR/OSCAR2008H.dat");\
  vector<cell*>* cell_array;

  // Initialize the scattering dynamics class
  scattering* dynamics = new scattering(&paraRdr);

  // Create a hard probe (hp)
  pdata hp_data;
  hp_data.particle_id = 1;
  hp_data.position += 0.6, 0.0, 0.0, 0.0;
  hp_data.velocity += 2.29416, 2.06474, 0.0, 0.0;
  hp_data.momentum += 0.688248, 0.619423, 0.0, 0.0;
  particle* hp = new particle(&paraRdr, &hp_data);

  // Time loop
  //========================================================
   for(int imin=0; imin<300; imin+=100){
    oscar->populateOSCAR(&paraRdr,imin,imin+100);
    cell_array = oscar->getHydro();
    for(int itime=imin; itime<(imin+100); itime+=1){

      // get hard probe (hp) coordinates and fluid cell (fc) container
      vector<int> hp_coord = hp->get_coord();
      int ifc = oscar->grabCellfromCoord(hp_coord); 
      cell* fc = (*cell_array)[ifc]; 
      cout << "coord: " << vector2string(hp_coord) << endl;

      // boost hard probe (hp) to the rest frame of the fluid cell (fc)
      vector<double> fc_vel = fc->velocity;
      cout << "fc_vel: " << vector2string(fc_vel) << endl;
      vector<double> p1_lab = hp->get_momentum();
      cout << "p1_lab: " << vector2string(p1_lab) << endl;
      vector<double> p1_cell = hp->lorentzboost(p1_lab,fc_vel);
      cout << "p1_cell: " << vector2string(p1_cell) << endl;
      hp->set_momentum(p1_cell);

      // sample a thermal quark in the fluid cell (fc) rest frame
      double mass = 0.2, temp = fc->thermal[2];
      //vector<double> p2_cell = oscar->sample_boltzmann(mass, temp);
      vector<double> p2_cell;
      p2_cell += 0.688248, -0.619423, 0.0, 0.0;
      cout << "p2_cell: " << vector2string(p2_cell) << endl;

      pdata tq_data;
      tq_data.particle_id = 1;
      hp_data.coord = hp_coord;
      tq_data.momentum = p2_cell;
      particle* tq = new particle(&paraRdr, &tq_data);

      // find the cms velocity of the collision pair
      vector<double> u_cms = oscar->getUcms(p1_cell, p2_cell);
      cout << "u_cms: " << vector2string(u_cms) << endl;

      // boost p1 and p2 to their cms frame
      vector<double> p1_cms = hp->lorentzboost(p1_cell, u_cms);
      cout << "p1_cms: " << vector2string(p1_cms) << endl;
      vector<double> p2_cms = tq->lorentzboost(p2_cell, u_cms);
      cout << "p2_cms: " << vector2string(p2_cms) << endl;
      hp->set_momentum(p1_cms);
      tq->set_momentum(p2_cms);
     
      // rotate p1 and p2 such that p1 lies along the z-axis
      double psi1 = hp->getPsi(p1_cms);
      double phi1 = hp->getPhi(p1_cms);
      vector<double> p1_cms_rot = hp->rotate(p1_cms,psi1,phi1);
      cout << "p1_cms_rot: " << vector2string(p1_cms_rot) << endl;
      vector<double> p2_cms_rot = tq->rotate(p2_cms,psi1,phi1);
      cout << "p2_cms_rot: " << vector2string(p2_cms_rot) << endl;
      hp->set_momentum(p1_cms_rot);
      tq->set_momentum(p2_cms_rot);

      // scatter p1 and p2, returning deflected theta and phi
      double s = dynamics->getMandelstamS(p1_cms_rot, p2_cms_rot);
      vector<int> ij2kl = dynamics->sample2to2(s,temp,1);
      double theta13 = dynamics->sampleTheta(s,temp,ij2kl);
      cout << "theta13: " << theta13 << endl;
      double phi13 = dynamics->samplePhi();
      cout << "phi13: " << phi13 << endl;

      // invert rotation
      //vector<double> p3_cms = hp->rotate(p1_cms,psi1,phi1);

      // increment quark evolution
      //hp->printPosition("data/quark_position_0.dat");
      //hp->stream(0.02);

      // increment hydro evolution
      //oscar->populate2D(itime,ithermal); 
      //oscar->print2D(path);  
      delete tq;
    }
    }
   //===============================================
  delete oscar;
  delete dynamics;
  return 1;
}
/*

  pdata quarkdata1;
  quarkdata1.particle_id += 1;
  quarkdata1.position += 0.6, -1, -2.5, 0.0;
  quarkdata1.velocity += 1.1547, 0.0, 0.57735, 0.0;
  particle* quark1 = new particle(&paraRdr, &quarkdata1);

  pdata quarkdata2;
  quarkdata2.particle_id += 1;
  quarkdata2.position += 0.6, 0, -1, 0.0;
  quarkdata2.velocity += 1.40028, 0.693103, 0.693103, 0.0;
  particle* quark2 = new particle(&paraRdr, &quarkdata2);

  pdata quarkdata3;
  quarkdata3.particle_id += 1;
  quarkdata3.position += 0.6, -2, 1, 0.0;
  quarkdata3.velocity += 1, 0.0, 0.0, 0.0;
  particle* quark3 = new particle(&paraRdr, &quarkdata3);

  pdata quarkdata4;
  quarkdata4.particle_id += 1;
  quarkdata4.position += 0.6, -1, -2.5, 0.0;
  quarkdata4.velocity += 1.1547, 0.0, -0.57735, 0.0;
  particle* quark4 = new particle(&paraRdr, &quarkdata4);

  quark1->stream(0.02);
  quark2->stream(0.02);
  quark3->stream(0.02);
  quark4->stream(0.02);

  quark1->printPosition("data/quark_position_1.dat");
  quark2->printPosition("data/quark_position_2.dat");
  quark3->printPosition("data/quark_position_3.dat");
  quark4->printPosition("data/quark_position_4.dat");

*/
