#include <iostream>

#include <sys/stat.h> // mkdir()
#include <cstring> // strerror()

#include "particle_data.h"
#include "peloton_solver.h"
#include "particle_viewer.h"

using namespace std;


int main(){
  //initialize the files

  
  // initialize the initializer  may use input file in later developement
  Initializer *init = new Initializer();

  // use initializer to initialize the global data
  Global_Data *gdata = new Global_Data(init);

  // initialize the viewer part, used to output data
  ParticleViewer *viewer = new ParticleViewer(gdata,outputfileNameFluid);

  // initialize the lp solver
  LPSolver *lpsolver = new LPSolver(init, gdata, viewer);

  // initialize the fluid particles, may add using restart later
  gdata->initFluidParticles_line();

  // solve the pde
  lpsolver->solver_1d();

  // clean the data
  gdata->cleanUpArrays();

  // clean the objects
  delete init;
  delete gdata;
  delete lpsolver;
  
  return 0;
}
