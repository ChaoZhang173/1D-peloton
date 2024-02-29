#include <iostream>
#include "lp_solver.h"

using namespace std;

LPSolver::LPSolver(Initializer *init, Global_Data *g, ParticleViewer *v) {
    gdata = g;
}

void LPSolver::solve_1d(){

    while(currenttime < tend)
    {
        // generate boundary particles

        // reorder particles, make sure particle position is in order

        // compute the cfl dt
        computeCFLCondition();
        // adjust dt by write time interval
        bool iswritestep = adjustDtByWriteTimeInterval(); 

        // heating model
        
        // lax-wendroff solver
        solve_laxwendroff();

        // compute boundary condition

        // update particle states

        // radiation cooling

        // move particles

        // coupute temperature

        // output data

        // delete remote particles

    }



}

void LPSolver::solve_laxwendroff(){
    const double *invelocity, *inpressure, *involume, *insoundspeed;
    double *outvelocity, *outpressure, *outvolume, *outsoundspeed;

    // spatial derivative: u: velocity, p: pressure, v: volume
    double ux, uxx, px, pxx, vx, vxx;
    
    // data for one particle
    pdata_t *pad;

    int li, lpnum = gdata->particle_data->elem_count;


  
}
