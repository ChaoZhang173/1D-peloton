#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include "initializer.h"
#include "boundary.h"

/*
This file is used for particle data. It has a struct pdata and a class Global_Data
pdata: a data structure for a single particle
Gloabal_Data::particle_data: a class that stores the list of all particle data (pdata) and functions
*/
typedef struct pdata{
  double x; //! coordinate
  double v; //! velocity
  double pressure;
  double soundspeed;
  double temperature;
  double volume;
  double sphdensity;
  double pressureT1; //! temperature used for updating states for each timestep
  double volumeT1;
  double soundspeedT1;
  double mass;

  bool ifboundary; //! if this is a boundary particle

  double leftintegral;
  double rightintegral;
  double rhointegral;
  double qplusminus;
} pdata_t; //! an alias for stuct pdata, usage: pdata_t p;

class Global_Data{
public:
  friend PelletSolver; //! PelletSolver could access private and proteced memebers
  Global_Data(Initializer* init);
  // destructor
  ~Global_Data();

  void initFluidParticles_line();

  Initializer *initializer;
  Geometry *geometry;
  State *state;
  EOS *eos;
  PelletSolver *pellet_solver;

  double gamma;
  double pinf;
  double einf;

  double backgroundpressure; 

  int pelletnumber; //! =1
  double currenttime;

  //--------print option-----------
  int printvelocity;
  int printdensity;
  int printmass;
  int printpressure;
  int printsoundspeed;
private:

};


#endif
