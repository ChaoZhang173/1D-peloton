#include <iostream>
#include <cmath>
#include <cassert>
#include "lp_solver.h"


using namespace std;

LPSolver::LPSolver(Initializer *init, Global_Data *g, ParticleViewer *v) {
    gdata = g;
    cfldt = 0.0;
    cflcoeff = init->getCFLcoeff();
    iswritestep = true;
    tstart = init->getTstart();
    tend = init->getTend();
    writetimeinterval = init->getWriteTimeInterval();
    currenttime = tstart;
    timestep = 0;
    nextwritetime = tstart + writetimeinterval;
    gdata->currenttime = currenttime;

    pellet_solver = new PelletSolver(init, gdata);
}

void LPSolver::solve_1d(){

    while(currenttime < tend)
    {
        // generate boundary particles
        for(int id = 0; id < gdata->boundarynumber; id++){
            gdata->boundary[id]->generateBoundaryParticle(gdata,gdata->eos,gdata->initialspacing,cfldt);
        }
        // reorder particles, make sure particle position is in order
        gdata->reorderParticles();
        // generate ghost partilces(for boundary particles)
        gdata->generateGhostParticles();
        // compute the cfl dt
        computeCFLCondition();
        // adjust dt by write time interval, update currenttime
        iswritestep = adjustDtByWriteTimeInterval(); 

        // heating model
        pellet_solver->heatingModel(currenttime);
        // lax-wendroff solver
        solve_laxwendroff();

        // compute boundary condition
        pellet_solver->computeBoundaryCondition(gdata, cfldt, gdata->initialspacing/5);
        // update particle states
        gdata->updateParticleStates();
        // update local spacing, need to be updated
        updateLocalSpacing();
        // radiation cooling

        // move particles
        moveParticle();       
        // coupute temperature
        computeTemperature();
        // output data

        // delete remote particles/ghost particles

        // update the ifStart
        gdata->ifStart = false;
    }



}

// ! the -∇·q term is implemented inside
void LPSolver::solve_laxwendroff() {
    const double *invelocity, *inpressure, *involume, *insoundspeed;
    double *outvelocity, *outpressure, *outvolume, *outsoundspeed;

    // spatial derivative: u: velocity, p: pressure, v: volume
    // order: /dx, /dx^2
    double Ud[2] = {0.,0.};
    double Pd[2] = {0.,0.};
    double Vd[2] = {0.,0.};
    
    // data for one particle
    pdata *pad;

    int li, lpnum = gdata->particle_data->size();

    for(li = 0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        
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
        computeSpatialDer(pad, Ud, Pd, Vd);

        // time integration
        timeIntegration(*invelocity, *inpressure, *involume, *insoundspeed, Ud, 
                                Pd, Vd, outvelocity, outpressure, outvolume);

        // check if time integration gives nan value
        if(isnan(*outvelocity) || isnan(*outpressure) || isnan(*outvolume)) {
        cout<<"[timeIntegration] Detect a particle has nan value!"<<endl;
        cout<<"[timeIntegration] Particle ID= "<<pad->id<<", position = "<<(pad->x)<<endl;
        assert(false);
        }

        // add -∇·q
        *outpressure += cfldt*pad->deltaq*((*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure)-1);

        // coumpute soundspeed
        *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume));

        // assign vaules to out state
        pad->oldv = *outvelocity;
        pad->pressureT1 = *outpressure;
        pad->volumeT1 = *outvolume;
        pad->soundspeedT1 = *outsoundspeed;
        }
}


void LPSolver::computeSpatialDer(pdata *pad, double *Ud, double *Pd, double *Vd) {

    pdata *pad_neil = &((*pad->neighbourparticle)[0]);
    pdata *pad_neir = &((*pad->neighbourparticle)[1]);

    // the coefficients for Newton Interpolation
    double coefP[6] ={0.,0.,0.,0.,0.,0.};
    double coefV[6] ={0.,0.,0.,0.,0.,0.};
    double coefU[6] ={0.,0.,0.,0.,0.,0.};
    // pass array to function = pass the pointer to the first element of the array
    computeDivdDifference(pad, coefU, coefP, coefV);

    // compute spatial derivative
    Ud[0] = coefU[3]+(pad->x - pad_neil->x)*coefU[5];
    Ud[1] = 2*coefU[5];
    Pd[0] = coefP[3]+(pad->x - pad_neil->x)*coefP[5];
    Pd[1] = 2*coefP[5];
    Vd[0] = coefV[3]+(pad->x - pad_neil->x)*coefV[5];
    Vd[1] = 2*coefV[5];

}


void LPSolver::computeDivdDifference(pdata *pad, double *coefU, double *coefP, double *coefV) {

    pdata *pad_neil = &((*pad->neighbourparticle)[0]);
    pdata *pad_neir = &((*pad->neighbourparticle)[1]); 

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

void LPSolver::timeIntegration(double inVelocity, double inPressure, double inVolume, 
                  double inSoundspeed, double *Ud, double *Pd, double *Vd, double *outVelocity, 
                  double *outPressure, double *outVolume) {

    double dt = cfldt;
    double gamma = inSoundspeed*inSoundspeed/inVolume/inPressure;
    double pinf = gdata->pinf;

    double ux = Ud[0];
    double uxx = Ud[1];
    double px = Pd[0];
    double pxx = Pd[1];
    double vx = Vd[0];
    double vxx = Vd[1];

    *outVelocity = inVelocity - dt*inVolume*px + 0.5*dt*dt*(gamma*inVolume*inPressure*uxx - inVolume*ux*px);
    *outPressure = inPressure - dt*gamma*inPressure*ux 
                + 0.5*dt*dt*(gamma*gamma/inPressure*ux*ux + gamma*vx*inPressure/px+gamma*inVolume*inPressure*pxx
                + gamma*inPressure*ux*ux);
    *outVolume = inVolume +dt*inVolume*ux + 0.5*dt*dt*(-inVolume*vx*px - inVolume*inVolume*pxx);
}

void LPSolver::moveParticle() {
    pdata *pad;
    int li, lpnum = gdata->particle_data->size();
    double dt = cfldt;

    for(li = 0; li < lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->x += 0.5*dt*(pad->v + pad->oldv);
    }
}

void LPSolver::computeTemperature() {
    pdata *pad;
    int li, lpnum = gdata->particle_data->size();

    for(li = 0; li < lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->temperature = gdata->eos->getTemperature(pad->pressure, 1./pad->volume);
    }
}

void LPSolver::updateLocalSpacing(){
    pdata *pad;
    int li, lpnum = gdata->particle_data->size();

    for(li=0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->localspacing *= pad->volume/pad->volumeT1;
    }
}

void LPSolver::computeCFLCondition(){
    pdata *pad;
    int li, lpnum = gdata->particle_data->size();

    double mindt = 100;
    double dt;
    double soundspeed;

    for (li = 0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        soundspeed = pad->soundspeed;
        dt = pad->localspacing/soundspeed;
        if (dt<mindt) mindt = dt;

        cfldt = mindt * cflcoeff;   
    }
}

bool LPSolver::adjustDtByWriteTimeInterval(){
    if(currenttime + cfldt >= nextwritetime){
        cfldt = nextwritetime - currenttime;
        currenttime = nextwritetime;
        nextwritetime += writetimeinterval;
        // warning message if cfldt <0
        if(cfldt < 0){
            cout<<"[LPSolver] Warning: cfldt < 0!"<<endl;
            cout<<"[LPSolver] currenttime = "<<currenttime<<", nextwritetime = "<<nextwritetime<<endl;
            cout<<"[LPSolver] cfldt = "<<cfldt<<endl;
            assert(false);
        }
        return true;
    }
    else {
        currenttime += cfldt;
        return false;
    }
}
