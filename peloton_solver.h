#ifndef __PELOTON_SOLVER_H__
#define __PELOTON_SOLVER_H__

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

  void computeSpatialDer(); 
  void timeIntegration();

  void computeCFLCondition();

  void computeTemperature();

  double tstart;
  double tend;
  double currenttime;
  double cfldt;

  int timestep;
  double writetimeinterval;
};


#endif
