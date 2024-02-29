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
    double Ud[2] = {0.,0.};
    double Pd[2] = {0.,0.};
    double Vd[2] = {0.,0.};
    
    // data for one particle
    pdata_t *pad;

    int li, lpnum = gdata->particle_data->size();

    for(li = 0; li<lpnum; li++){
        pad = &(gdata->particle_data[li]);
        
        // if the particle is at the boundray (inside the pellet), skip it
        if (pad->ifboundary) 
            continue;
        
        // get the in state
        invelocity = &(pad->v);
        inpressure = &(pad->pressure);
        involume = &(pad->volume);
        insoundspeed = &(pad->soundspeed);

        // set out state 
        outvelocity = &(pad->oldv);
        outpressure = &(pad->pressureT1);
        outvolume = &(pad->volumeT1);
        outsoundspeed = &(pad->soundspeedT1);

        // set out state = in state (incase solver fails)
        *outvelocity = *invelocity;
        *outpressure = *inpressure;
        *outvolume = *involume;
        *outsoundspeed = *insoundspeed;

        // if the particle has 0 volume or 0 soundspeed, skip it
        if(*insoundspeed < 1e-10 || *involume < 1e-10) {
            cout<<"[LW] Detect a particle has 0 volume or 0 soundspeed!"<<endl;
            cout<<"[LW] Particle ID= "<<li<<", position = "<<(pad->x)<<endl;
            continue;
        }

        // compute spatial derivative
        computeSpatialDer(pad, invelocity, inpressure, involume, Ud, Pd, Vd);

        // time integration
        timeIntegration();

        // assign vaules to out state
    }

  
}
