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
    // order: /dx, /dx^2
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

void LPSolver::computeSpatialDer(pdata_t *pad, double *Ud, double *Pd, double *Vd) {
    // next: add null pointer check
    pdata_t *pad_neil = &((*pad->neighbourparticle)[0]);
    pdata_t *pad_neir = &((*pad->neighbourparticle)[1]);

    // the coefficients for Newton Interpolation
    double *coefP[6] =[0.,0.,0.,0.,0.,0.];
    double *coefV[6] =[0.,0.,0.,0.,0.,0.];
    double *coefU[6] =[0.,0.,0.,0.,0.,0.];

    computeDivdDifference(pad, coefU, coefP, coefV);

    // compute spatial derivative
    Ud[0] = coefU[3]+(2*pad->x - pad_neil->x - pad_neir->x)*coefU[5];
    Ud[1] = 2*coefU[5];
    Pd[0] = coefP[3]+(2*pad->x - pad_neil->x - pad_neir->x)*coefP[5];
    Pd[1] = 2*coefP[5];
    Vd[0] = coefV[3]+(2*pad->x - pad_neil->x - pad_neir->x)*coefV[5];
    Vd[1] = 2*coefV[5];

}


void LPSolver::computeDivdDifference(pdata_t *pad, double *coefU, double *coefP, double *coefV) {
    // next: add null pointer check
    pdata_t *pad_neil = &((*pad->neighbourparticle)[0]);
    pdata_t *pad_neir = &((*pad->neighbourparticle)[1]); 

    coefU[0] = pad_neil->v;
    coefU[1] = pad->v;
    coefU[2] = pad_neir->v;
    coefU[3] = (coefU[1] - coefU[0])/(pad->x - pad_neil->x);
    coefU[4] = (coefU[2] - coefU[1])/(pad_neir->x - pad->x);
    coefU[5] = (coefU[4] - coefU[3])/(pad_neir->x - pad_neil->x);

    coefP[0] = pad_neil->pressure;
    coefP[1] = pad->pressure;
    coefP[2] = pad_neir->pressure;
    coefP[3] = (coefP[1] - coefP[0])/(pad->x - pad_neil->x);
    coefP[4] = (coefP[2] - coefP[1])/(pad_neir->x - pad->x);
    coefP[5] = (coefP[4] - coefP[3])/(pad_neir->x - pad_neil->x);

    coefV[0] = pad_neil->volume;
    coefV[1] = pad->volume;
    coefV[2] = pad_neir->volume;
    coefV[3] = (coefV[1] - coefV[0])/(pad->x - pad_neil->x);
    coefV[4] = (coefV[2] - coefV[1])/(pad_neir->x - pad->x);
    coefV[5] = (coefV[4] - coefV[3])/(pad_neir->x - pad_neil->x);
    
}