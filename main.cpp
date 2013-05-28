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
  pdata p1_data, p2_data, p3_data, p4_data;
  vector<int> p1_coord, ij2kl;
  vector<double> fc_vel, p1_lab, p1_cell, p2_cell, u_cms, p1_cms, p2_cms;
  vector<double> p3_cms, p4_cms, u_cms_inv, p3_lab, p4_lab;
  double mass, temp, s, theta13, phi13, E3, E4;
  int ifc;
  cell* fc;

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

  // Create a hard probe (p1)
  p1_data;
  p1_data.particle_id = 1;
  p1_data.position += 0.6, 0.0, 0.0, 0.0;
  p1_data.velocity += 2.29416, 2.06474, 0.0, 0.0;
  p1_data.momentum += 0.688248, 0.619423, 0.0, 0.0;
  particle* p1 = new particle(&paraRdr, &p1_data);

  // Time loop
  //========================================================
   for(int imin=0; imin<300; imin+=20){
    oscar->populateOSCAR(&paraRdr,imin,imin+20);
    cell_array = oscar->getHydro();
    for(int itime=imin; itime<(imin+20); itime+=1){

      //###############################################################
      //############## Execute 2->2 Scattering Event ##################

      // get hard probe (p1) coordinates and fluid cell (fc) container
      //==============================================================
      p1_coord = p1->get_coord();
      ifc = oscar->grabCellfromCoord(p1_coord); 
      fc = (*cell_array)[ifc]; 
      
      // boost hard probe (p1) to the rest frame of the fluid cell (fc)
      //==============================================================
      fc_vel = fc->velocity;
      p1_lab = p1->get_momentum();
      p1_cell = p1->lorentzboost(p1_lab,fc_vel);
      p1->set_momentum(p1_cell);

      // sample a thermal quark (p2) in the fluid cell (fc) rest frame
      //=============================================================
      //vector<double> p2_cell = oscar->sample_boltzmann(mass, temp);
      mass = 0.2;
      temp = fc->thermal[2];
      p2_cell += 0.688248, -0.619423, 0.0, 0.0;
      p2_data.particle_id = 1;
      p2_data.coord = p1_coord;
      p2_data.momentum = p2_cell;
      particle* p2 = new particle(&paraRdr, &p2_data); 

      // boost p1 and p2 to their cms frame
      //============================================================
      u_cms = oscar->getUcms(p1_cell, p2_cell);
      p1_cms = p1->lorentzboost(p1_cell, u_cms);
      p2_cms = p2->lorentzboost(p2_cell, u_cms);
      p1->set_momentum(p1_cms);
      p2->set_momentum(p2_cms);

      // scatter p1 and p2, returning deflected theta13 and phi13
      //===========================================================
      s = dynamics->getMandelstamS(p1_cms, p2_cms);
      ij2kl = dynamics->sample2to2(s,temp,1);
      theta13 = dynamics->sampleTheta(s,temp,ij2kl);
      phi13 = dynamics->samplePhi();

      // create deflected particles p3 and p4
      //==========================================================
      p3_data.particle_id = ij2kl[2]; 
      p4_data.particle_id = ij2kl[3];
      p3_data.coord = p1->get_coord();
      p4_data.coord = p1->get_coord();
      particle* p3 = new particle(&paraRdr, &p3_data);
      particle* p4 = new particle(&paraRdr, &p4_data);
      p3_cms = p1->rotate(p1_cms, M_PI-theta13, phi13);
      p4_cms = p2->rotate(p2_cms, theta13-M_PI, phi13);
      p3->set_momentum(p3_cms);
      p4->set_momentum(p4_cms);

      // restore p3 and p4 to the lab frame
      //=========================================================
      u_cms_inv = p3->reflectfourvector(u_cms);
      p3_lab = p3->lorentzboost(p3_cms, u_cms_inv);
      p4_lab = p4->lorentzboost(p4_cms, u_cms_inv);
      p3->set_momentum(p3_lab);
      p4->set_momentum(p4_lab);

      // follow the outgoing particle with the leading energy
      //========================================================
      E3 = p3_lab[0];
      E4 = p4_lab[0];
      if(E3>E4) p1 = p3;
      else if(E4>E3) p1 = p4;
      else if(drand48()>0.5) p1 = p3;
      else p1 = p4;

      delete p2; 
      delete p3; 
      delete p4;

      // free stream leading partons until next scattering event
      //=========================================================
      p1->printPosition("data/quark_position_0.dat");
      p1->stream(0.02);

      cout << "p1 position: " << vector2string(p1->get_position()) << endl;
      cout << "p1 momentum: " << vector2string(p1->get_momentum()) << endl;
      cout << endl;
      cout << endl;

      // increment hydro evolution
      //=========================================================
      //oscar->populate2D(itime,ithermal); 
      //oscar->print2D(path);  
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
