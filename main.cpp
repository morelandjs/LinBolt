#include <iostream>
#include <stdio.h>
#include "system.h"
#include "OSCAR.h"
#include "particle.h"
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

  // Create a quark
  pdata quarkdata0;
  quarkdata0.particle_id += 1;
  quarkdata0.position += 0.6, 0.0, 0.0, 0.0;
  //quarkdata0.velocity += 1.0, 0.0, 0.0, 0.0;
  quarkdata0.velocity += 2.29416, 2.06474, 0.0, 0.0;
  particle* quark0 = new particle(&paraRdr, &quarkdata0);

  // Time loop
  //========================================================
   for(int imin=0; imin<300; imin+=100){
    oscar->populateOSCAR(&paraRdr,imin,imin+100);
    cell_array = oscar->getHydro();
    for(int itime=imin; itime<(imin+100); itime+=1){
      int label = itime;
      //sprintf(path,"data/energy_dens_%05d.dat",label);
      vector<int> quark_coord = quark0->get_coord();
      cout << "quark coord: " << vector2string(quark_coord) << endl;
      int cell_number = oscar->grabCellfromCoord(quark_coord); 
      cout << "cell number: " << cell_number << endl;
      cout << "size of cell array: " << (*cell_array).size() << endl;
      cell* local_cell = (*cell_array)[cell_number]; 
      vector<double> cell_velocity = local_cell->velocity;
      cout << "cell velocity: " << vector2string(cell_velocity) << endl;

      double mass = 0.0;
      double temperature = local_cell->thermal[2];
      cout << "temperature: " << temperature << endl;
      vector<double> momentum0 = oscar->sample_boltzmann(mass, temperature);
      cout << "cell momentum0: " << vector2string(momentum0) << endl;

      vector<double> boost;
      boost += cell_velocity[0], cell_velocity[1], cell_velocity[2], cell_velocity[3];
      cout << "boost: " << vector2string(boost) << endl;
      vector<double> momentum = oscar->lorentzboost(momentum0,boost);
      cout << "check boost: " << vector2string(oscar->lorentzboost(cell_velocity,boost)) << endl;
      cout << " cell momentum: " << vector2string(momentum) << endl;
 
      // increment quark evolution
      
      //quark0->printPosition("data/quark_position_0.dat");
     
      //free stream
      quark0->stream(0.02);

      // increment hydro evolution
      //oscar->populate2D(itime,ithermal); 
      //oscar->print2D(path);  
    }
    }
  //=======================================================

  // Delete the hydro array
  delete oscar;
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
