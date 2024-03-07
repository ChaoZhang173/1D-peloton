#ifndef __PELLET_SOLVER_H__
#define __PELLET_SOLVER_H__

#include <vector>

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

  void heatingModel(double currenttime);

  //! compute electron density integral using trapezoid rule, 
  // also set -∇·q and q+- to 0 for following computeHeatDeposition
  void computeDensityIntegral();
  //! with multiple heating source
  void computeHeatDeposition(double currenttime);
  //! add heating source
  void addHeatSource(double teinf, double neinf, double currenttime);
  //！ set parameters of pellet material
  void setPelletMaterial(int materialid);

  int heatsource_numer;//! the number of heating sources
  std::vector<double> teinf;
  std::vector<double> neinf;
  double warmuptime;

  double mu;
  double mass;
  double Z;
  double I;
  double sublimationenergy;
  double neeff;

  double one_plus_Zstar; // ! 1 + Z*

  Global_Data *gdata;

  const double heatK = 1.602e-18; // ! 1.602e-18, 1eV = 1.602e-19J
  const double masse = 9.109e-28; // ! mass of electron, 9.109e-28g
};
#endif
