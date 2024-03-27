#ifndef __PARTICLE_DATA_H__
#define __PARTICLE_DATA_H__

#include <vector>
#include <memory>

#include "initializer.h"
#include "boundary.h"
#include "state.h"

/*
This file is used for particle data. It has a struct pdata and a class Global_Data
pdata: a data structure for a single particle
Gloabal_Data::particle_data: a class that stores the list of all particle data (pdata) and functions
*/
struct pdata{
  double x; //! coordinate
  double v; //! velocity
  double oldv; //! velocity used to update states for each timestep
  double pressure;
  double soundspeed;
  double temperature;
  double volume;
  //double sphdensity;// currently not used
  double pressureT1; //! temperature used for updating states for each timestep
  double volumeT1;
  double soundspeedT1;
  double mass;

  // need to be initialized
  double localspacing; //! 0.5(dis_left + dis_right)

  bool ifboundary; //! if this is a boundary particle
  
  //! for 1D code, only has left integral,
  // it is from the right neighbour to the particle itself
  double leftintegral; 
  double rightintegral; // ! not used
  double rhointegral;
  double qplusminus;
  double deltaq; // ! -∇·q

  //! neighbor particles
  pdata *leftneighbour = nullptr;
  pdata *rightneighbour = nullptr;

  // neon radiation cooling
  double radcool;

};

class PelletSolver;

class Global_Data{
public:
  friend class PelletSolver; //! PelletSolver could access private and proteced memebers
  Global_Data(Initializer* init);
  // destructor
  ~Global_Data();


  //! initialize fluid particles and ghost particles 
  // the particles are initialized in a layer that close to the pellet surface 
  // the first particle is set as the pellet surface,
  // all the particles are set as fluid particles
  void initFluidParticles_line();

  //! update particle states for each timestep
  void updateParticleStates();

  //! set the eos, eoschoice = 1 for polytropic, 2 for saha neon 
  void setEOS();

  //! reorder particles and set neighbours
  void reorderParticles();

  //! generate ghost particles
  void generateGhostParticles();

  //! update local spacing: 0.5((distance to left) + (distance to right))
  // if two particles too close, replace the particle
  void updatelocalSpacing();

  //! delete particles with bad states(small pressure, large volume, high temperature, etc.)
  void deleteBadParticles();

  //! true: at the first step
  bool ifStart;

  Initializer *initializer;
  State *state;
  EOS *eos;
  PelletSolver *pellet_solver;

  //! 1: polytropic; 2: saha neon; 
  int eoschoice;

  double gamma;
  double pinf;
  double einf;

  double initialspacing;
  
  // smallest dx between particles
  double mindx; 
  
  double backgroundpressure; 

  // invalid pressure and density
  double invalidpressure;
  double invaliddensity;

  // bad pressure and volume, used to delete bad particles
  double badpressure;
  double badvolume;

  int pelletnumber; //! =1
  double currenttime;

  //! the length of layer that will generate initial particles
  double initiallayerlength;

  // a vector that stores all the particle data
  std::unique_ptr<std::vector<pdata>> particle_data;
  // a vector that stores the ghost particles
  std::unique_ptr<std::vector<pdata>> ghostparticle_data;

  //! the estimated max number of particles, 
  // used to preallocate memory for particle list
  // could be changed if the size is not enough
  // currently using smart pointer, not used
  int maxparticlenum;

  //! boundary objects
  std::vector<Boundary*> boundary;
  //! number of boundary objects
  int boundarynumber;

  //--------print option-----------
  int printvelocity;
  int printdensity;
  int printmass;
  int printpressure;
  int printsoundspeed;
private:

};


#endif
