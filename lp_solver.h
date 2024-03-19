#ifndef __LP_SOLVER_H__
#define __LP_SOLVER_H__

#include "initializer.h"
#include "particle_data.h"
#include "particle_viewer.h"
#include "pellet_solver.h"


/*
This file is for LPSolver class, it implements functions to solve PDE using lax-wendroff method
*/
class LPSolver {
public:
  LPSolver(Initializer* init, Global_Data *g, ParticleViewer *viewer);

  ~LPSolver(){}
  
  void solve_1d();

  void moveParticle();
  // for all particles expect right boundary
  void solve_laxwendroff();
  // for right boundary, using upwind
  void solve_upwind_right_boundary();
  //! compute spatial derivatives using Newton Interpolation
  // order: /dx, /dx^2
  void computeSpatialDer(pdata *pad, double *Ud, double *Pd, double *Vd); 
  //! compute spatial derivatives for upwind,
  // order: /dx_left, /dx_right
  void computeSpatialDer_upwind(pdata *pad, double *Ud, double *Pd, double *Vd);
  //! compute spatial derivatives using devided difference
  void computeSpatialDer_DD(pdata *pad, double *Ud, double *Pd, double *Vd);
  //! time integration using Lax-Wendroff method
  void timeIntegration(double inVelocity, double inPressure, double inVolume, 
                  double inSoundspeed, double *Ud, double *Pd, double *Vd, double *outVelocity, 
                  double *outPressure, double *outVolume);
  //! time integration using upwind method
  void timeIntegration_upwind(double inVelocity, double inPressure, double inVolume, 
                  double inSoundspeed, double *Ud, double *Pd, double *Vd, double *outVelocity, 
                  double *outPressure, double *outVolume);              

  void computeCFLCondition();

  void computeTemperature();
  //! update localspacing using volumeT1, used in 3D code, not used in 1D code
  void updateLocalSpacing();
  // ! update currenttime, nextwritetime and gdata->currenttime 
  bool adjustDtByWriteTimeInterval();

  
  //! compute the coefficients for Newton Interpolation, using divided difference
  //  coefF = f[t1], f[t2], f[t3], f[t1,t2], f[t2,t3], f[t1,t2,t3]
  void computeDivdDifference(pdata *pad, double *coefU, double *coefP, double *coefV);

  Global_Data *gdata; 
  ParticleViewer *viewer;
  PelletSolver *pellet_solver;

  double tstart;
  double tend;
  double currenttime;
  double cfldt;
  double cflcoeff;
  double nextwritetime;

  int timestep;
  double writetimeinterval;

  bool iswritestep;
  double writestep;
};


#endif
