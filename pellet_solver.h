#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__

#include "particle_data.h"
#include "material_lib.h"
#include "initializer.h"

/*
This file is for pellet solver.
The class PelletSolver stores all functions that related to pellet problem.
*/

class PelletSolver{
public:
  PelletSolver (Initializer* init,Global_Data *gdata);
  ~PelletSolver();

  void heatingModel(double dt);


  double mu;
  double mass;
  double Z;
  double I;
  double sublimationenergy;
  double teinf;
  double neinf;
  double neeff;
};
#endif
