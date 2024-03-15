#include <iostream>

#include <sys/stat.h> // mkdir()
#include <cstring> // strerror()

#include "particle_data.h"
#include "lp_solver.h"
#include "particle_viewer.h"

using namespace std;


int main(int argc, char *argv[]){
  //initialize the files
  if (argc != 2){
    cerr <<"[Error] Invalid number of input arguments.\n";
    cerr << "Usage: " << argv[0] << " <output file name>" << endl;
    return 1;
  }
  string outputFilename = argv[1];
  // initialize the initializer  may use input file in later developement
  Initializer *init = new Initializer();

  // use initializer to initialize the global data
  Global_Data *gdata = new Global_Data(init);

  // initialize the viewer part, used to output data
  ParticleViewer *viewer = new ParticleViewer(gdata,outputFilename);

  // initialize the lp solver
  LPSolver *lpsolver = new LPSolver(init, gdata, viewer);

  // initialize the fluid particles, and ghost particles may add using restart later
  gdata->initFluidParticles_line();

  // solve the pde
  lpsolver->solve_1d();

  // clean the data
  

  // clean the objects
  delete init;
  delete gdata;
  delete lpsolver;
  delete viewer;
  
  return 0;
}
