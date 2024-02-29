#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__

#include "initializer.h"
#include "particle_data.h"
#include "particle_viewer.h"
#include "pellet_solver.h"

/*
This file is for LPSolver class, it implements functions to solve PDE using lax-wendroff method
*/
clas LPSolver {
public:
  LPSolver(Initializer* init, Global_Data *g, ParticleViewer *viewer);

  ~LPSolver(){}
  
  void solve_1d();

  void moveParticle();
  void solve_laxwendroff();

  void computeSpatialDer(pdata_t *pad, const double *inVelocity, const double *inPressure, 
  const double *inVolume, double *Ud, double *Pd, double *Vd); 
  void timeIntegration();

  void computeCFLCondition();

  void computeTemperature();

  Global_Data *gdata; 
  ParticleViewer *viewer;
  PelletSolver *pellet_solver;

  double tstart;
  double tend;
  double currenttime;
  double cfldt;

  int timestep;
  double writetimeinterval;
};


#endif
