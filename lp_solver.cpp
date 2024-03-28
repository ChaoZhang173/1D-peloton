#include <iostream>
#include <cmath>
#include <cassert>
#include "lp_solver.h"


using namespace std;

LPSolver::LPSolver(Initializer *init, Global_Data *g, ParticleViewer *v) {
    gdata = g;
    viewer = v; 
    cfldt = 0.0;
    cflcoeff = init->getCFLcoeff();
    iswritestep = true;
    writestep = -1;
    tstart = init->getTstart();
    tend = init->getTend();
    writetimeinterval = init->getWriteTimeInterval();
    currenttime = tstart;
    timestep = 0;
    nextwritetime = tstart;
    gdata->currenttime = currenttime;

    pellet_solver = new PelletSolver(init, gdata);
}

void LPSolver::solve_1d(){

    while(currenttime <= tend)
    {
        // generate boundary particles
        for(int id = 0; id < gdata->boundarynumber; id++){
            gdata->boundary[id]->generateBoundaryParticle(gdata,gdata->eos,gdata->initialspacing,cfldt,&(pellet_solver->mass_accum));
        }
        // reorder particles, make sure particle position is in order, set neighbours
        gdata->reorderParticles();
        // generate(update) ghost partilces(for boundary particles)
        gdata->generateGhostParticles();
        // update loacal spacing, replace particles if too close
        gdata->updatelocalSpacing();
        // compute the cfl dt
        computeCFLCondition();
        // adjust dt by write time interval, update currenttime
        iswritestep = adjustDtByWriteTimeInterval(); 
        if(iswritestep) {
            cout<<"[LPSolver] currenttime = "<<currenttime<<", nextwritetime = "<<nextwritetime<<endl;
            viewer->writeResult(writestep, currenttime);
        }
        // heating model
        pellet_solver->heatingModel(currenttime);
//        // lax-wendroff solver for other particles
//        solve_laxwendroff();
//        // calculate right & left boundaries using upwind method
//        solve_upwind_boundries();

        // only use upwind method for now
        solve_upwind();

        // compute boundary condition, the choice of dx needs to be updated
        pellet_solver->computeBoundaryCondition(gdata, cfldt, gdata->initialspacing);
        // update particle states
        gdata->updateParticleStates();
        // radiation cooling
        pellet_solver->neonRadiationCooling(cfldt);
        // move particles
        moveParticle();  
        // update local spacing based on volume change    
        //updateLocalSpacing();
        // reorder particle and update loacl spacing, move to next timestep
        //gdata->updatelocalSpacing();
        // coupute temperature
        computeTemperature();
        // output data, move to begining of the next timestep

        // delete remote particles/ghost particles
        gdata->deleteBadParticles();
        // update the ifStart
        gdata->ifStart = false;
    }



}

// ! the -∇·q term is implemented inside
void LPSolver::solve_laxwendroff() {
    cout<<"[LPSolver] Entering solve_laxwendroff..."<<endl;
    const double *invelocity, *inpressure, *involume, *insoundspeed;
    double *outvelocity, *outpressure, *outvolume, *outsoundspeed;

    // spatial derivative: u: velocity, p: pressure, v: volume
    // order: /dx, /dx^2
    double Ud[2] = {0.,0.};
    double Pd[2] = {0.,0.};
    double Vd[2] = {0.,0.};
    
    // data for one particle
    pdata *pad;

    size_t li, lpnum = gdata->particle_data->size();
    // for all particles except the right & left boundaries
    for(li = 1; li<lpnum-1; li++){
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
        //computeSpatialDer_DD(pad, Ud, Pd, Vd);

        // time integration
        timeIntegration(*invelocity, *inpressure, *involume, *insoundspeed, Ud, 
                                Pd, Vd, outvelocity, outpressure, outvolume);
        
        // add -∇·q
        *outpressure += cfldt*pad->deltaq*((*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure)-1);
        
        // check if time integration gives nan value
        if(isnan(*outvelocity) || isnan(*outpressure) || isnan(*outvolume)) {
        cout<<"[timeIntegration] Detect a particle has nan value!"<<endl;
        cout<<"[timeIntegration] Particle ID= "<<li<<", position = "<<(pad->x)<<endl;
        assert(false);
        }
        if(isinf(*outvelocity) || isinf(*outpressure) || isinf(*outvolume)) {
        cout<<"[timeIntegration] Detect a particle has inf value!"<<endl;
        cout<<"[timeIntegration] Particle ID= "<<li<<", position = "<<(pad->x)<<endl;
        assert(false);
        }
        // check if the value is within the physical range
        if(*outpressure < gdata->invalidpressure || 1./(*outvolume) < gdata->invaliddensity) {
            double gamma = (*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure);
            cout<<"[LW Solver] inpressure: "<<*inpressure<<", outpressure: "<<*outpressure<<endl;
            cout<<"   involume: "<<*involume<<", outvolume: "<<*outvolume<<endl;
            cout<<"   position: "<<pad->x<<"li: "<<li<<", heating: "<<(gamma-1)*cfldt*pad->deltaq<<endl;
            assert(false);
        }

        // coumpute soundspeed
        *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume), li);

        // assign vaules to out state
        pad->oldv = *outvelocity;
        pad->pressureT1 = *outpressure;
        pad->volumeT1 = *outvolume;
        pad->soundspeedT1 = *outsoundspeed;
        }
    cout<<"[LPSolver] solve_laxwendroff finished!"<<endl;
}


void LPSolver::computeSpatialDer(pdata *pad, double *Ud, double *Pd, double *Vd) {

    pdata *pad_neil = pad->leftneighbour;
    //pdata *pad_neir = pad->rightneighbour;

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

void LPSolver::computeSpatialDer_DD(pdata *pad, double *Ud, double *Pd, double *Vd) {

    pdata *pad_neil = pad->leftneighbour;
    pdata *pad_neir = pad->rightneighbour;

    // the coefficients for Newton Interpolation
    double coefP[6] ={0.,0.,0.,0.,0.,0.};
    double coefV[6] ={0.,0.,0.,0.,0.,0.};
    double coefU[6] ={0.,0.,0.,0.,0.,0.};
    // pass array to function = pass the pointer to the first element of the array
    computeDivdDifference(pad, coefU, coefP, coefV);

    // compute spatial derivative
    Ud[0] = 0.5 * (coefU[4] + coefU[3]);
    Ud[1] = 2.0 * (coefU[4] - coefU[3])/(pad_neir->x - pad_neil->x);
    Pd[0] = 0.5 * (coefP[4] + coefP[3]);
    Pd[1] = 2.0 * (coefP[4] - coefP[3])/(pad_neir->x - pad_neil->x);
    Vd[0] = 0.5 * (coefV[4] + coefV[3]);
    Vd[1] = 2.0 * (coefV[4] - coefV[3])/(pad_neir->x - pad_neil->x);
}

void LPSolver::computeDivdDifference(pdata *pad, double *coefU, double *coefP, double *coefV) {

    pdata *pad_neil = pad->leftneighbour;
    pdata *pad_neir = pad->rightneighbour; 

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
    //double pinf = gdata->pinf; // currently not used

    double ux = Ud[0];
    double uxx = Ud[1];
    double px = Pd[0];
    double pxx = Pd[1];
    double vx = Vd[0];
    //double vxx = Vd[1];

    *outVelocity = inVelocity - dt*inVolume*px + 0.5*dt*dt*(gamma*inVolume*inPressure*uxx - inVolume*ux*px);
    *outPressure = inPressure - dt*gamma*inPressure*ux + 0.5*dt*dt*(gamma*gamma*inPressure*ux*ux 
                + gamma*vx*inPressure*px+gamma*inVolume*inPressure*pxx + gamma*inPressure*ux*ux);
    *outVolume = inVolume +dt*inVolume*ux + 0.5*dt*dt*(-inVolume*vx*px - inVolume*inVolume*pxx);
}

void LPSolver::solve_upwind(){
    
    cout<<"[LPSolver] Entering solve_upwind..."<<endl;
    pdata *pad;
    size_t li, lpnum = gdata->particle_data->size();

    for(li = 0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        solve_upwind_particle(pad, li);
    }
    cout<<"[LPSolver] solve_upwind finished!"<<endl;

}

void LPSolver::solve_upwind_boundaries(){
    
    cout<<"[LPSolver] Entering solve_upwind..."<<endl;
    pdata *pad;
    double lpnum = gdata->particle_data->size();
    // solve left boundary
    pad = &((*gdata->particle_data)[0]);
    solve_upwind_particle(pad, 0);
    // solve right boundary
    pad = &((*gdata->particle_data)[lpnum-1]);
    solve_upwind_particle(pad, lpnum-1);
    cout<<"[LPSolver] solve_upwind finished!"<<endl;
}

void LPSolver::solve_upwind_particle(pdata *pad, int pnum) {

    const double *invelocity, *inpressure, *involume, *insoundspeed;
    double *outvelocity, *outpressure, *outvolume, *outsoundspeed;

    // spatial derivative: u: velocity, p: pressure, v: volume
    // order: /dx_left, /dx_right
    double Ud[2] = {0.,0.};
    double Pd[2] = {0.,0.};
    double Vd[2] = {0.,0.};

    //double dt = cfldt;

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
        cout<<"[UW] Detect a particle has 0 volume or 0 soundspeed!"<<endl;
        cout<<"[UW] Particle ID= "<<pnum<<", position = "<<(pad->x)<<endl;
        return;
    }

    // compute spatial derivative
    computeSpatialDer_upwind(pad, Ud, Pd, Vd);

    // time integration
    timeIntegration_upwind(*invelocity, *inpressure, *involume, *insoundspeed, Ud, 
                            Pd, Vd, outvelocity, outpressure, outvolume, pad->x);

    // add -∇·q
    *outpressure += cfldt*pad->deltaq*((*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure)-1);

    // check if time integration gives nan/inf value
    if(isnan(*outvelocity) || isnan(*outpressure) || isnan(*outvolume)) {
    cout<<"[timeIntegration] Detect a particle has nan value!"<<endl;
    cout<<"[timeIntegration] Particle ID= "<<pnum<<", position = "<<(pad->x)<<endl;
    cout<<"invelocity: "<<*invelocity<<", outvelocity: "<<*outvelocity<<endl;
    cout<<"inpressure: "<<*inpressure<<", outpressure: "<<*outpressure<<endl;
    cout<<"involume: "<<*involume<<", outvolume: "<<*outvolume<<endl;
    assert(false);
    }
    if(isinf(*outvelocity) || isinf(*outpressure) || isinf(*outvolume)) {
    cout<<"[timeIntegration] Detect a particle has inf value!"<<endl;
    cout<<"[timeIntegration] Particle ID= "<<pnum<<", position = "<<(pad->x)<<endl;
    cout<<"invelocity: "<<*invelocity<<", outvelocity: "<<*outvelocity<<endl;
    cout<<"inpressure: "<<*inpressure<<", outpressure: "<<*outpressure<<endl;
    cout<<"involume: "<<*involume<<", outvolume: "<<*outvolume<<endl;
    assert(false);
    }
    // check if the value is within the physical range
    if(*outpressure < gdata->invalidpressure || 1./(*outvolume) < gdata->invaliddensity) {
        double gamma = (*insoundspeed)*(*insoundspeed)/(*involume)/(*inpressure);
        cout<<"[UW Solver] inpressure: "<<*inpressure<<", outpressure: "<<*outpressure<<endl;
        cout<<"   involume: "<<*involume<<", outvolume: "<<*outvolume<<endl;
        cout<<"   position: "<<pad->x<<", heating: "<<(gamma-1)*cfldt*pad->deltaq<<endl;
        assert(false);
    }

    // coumpute soundspeed
    *outsoundspeed = gdata->eos->getSoundSpeed(*outpressure, 1./(*outvolume),pnum);

    // assign vaules to out state
    pad->oldv = *outvelocity;
    pad->pressureT1 = *outpressure;
    pad->volumeT1 = *outvolume;
    pad->soundspeedT1 = *outsoundspeed;

}

void LPSolver::computeSpatialDer_upwind(pdata *pad, double *Ud, double *Pd, double *Vd) {

    //pdata *pad_neil = pad->leftneighbour;
    //pdata *pad_neir = pad->rightneighbour;

    // the coefficients for Newton Interpolation
    double coefP[6] ={0.,0.,0.,0.,0.,0.};
    double coefV[6] ={0.,0.,0.,0.,0.,0.};
    double coefU[6] ={0.,0.,0.,0.,0.,0.};
    // pass array to function = pass the pointer to the first element of the array
    computeDivdDifference(pad, coefU, coefP, coefV);

    Ud[0] = coefU[3];
    Ud[1] = coefU[4];
    Pd[0] = coefP[3];
    Pd[1] = coefP[4];
    Vd[0] = coefV[3];
    Vd[1] = coefV[4];

}

void LPSolver::timeIntegration_upwind(double inVelocity, double inPressure, double inVolume, 
                  double inSoundspeed, double *Ud, double *Pd, double *Vd, double *outVelocity, 
                  double *outPressure, double *outVolume, double x) {
    
    double vt,pt,ut;
    double dt = cfldt;
    //double gamma = inSoundspeed*inSoundspeed/inVolume/inPressure;
    double K = inSoundspeed*inSoundspeed/inVolume/inVolume;
    //double pinf = gdata->pinf; // currently not used

    double uxl = Ud[0];
    double uxr = Ud[1];
    double pxl = Pd[0];
    double pxr = Pd[1];
    //double vxl = Vd[0];
    //double vxr = Vd[1];

    vt = 0.5*inVolume*(uxr+uxl) - 0.5*inVolume*sqrt(x/K)*(pxr-pxl) + 2/x*inVolume*inVelocity;
    ut = 0.5 * inVolume * sqrt(x * K) * (uxr - uxl) - 0.5 * x * inVolume * (pxr + pxl);
    pt = -0.5*inVolume*K*(uxr+uxl) + 0.5*inVolume*sqrt(x*K)*(pxr-pxl) - 2/x*inVolume*K*inVelocity;

    *outVelocity = inVelocity + dt*ut;
    *outPressure = inPressure + dt*pt;
    *outVolume = inVolume + dt*vt;
}

void LPSolver::moveParticle() {
    pdata *pad;
    size_t li, lpnum = gdata->particle_data->size();
    double dt = cfldt;

    for(li = 0; li < lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->x += 0.5*dt*(pad->v + pad->oldv);
    }
}

void LPSolver::computeTemperature() {
    cout<<"[LPSolver] Entering computeTemperature..."<<endl;
    pdata *pad;
    size_t li, lpnum = gdata->particle_data->size();

    for(li = 0; li < lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->temperature = gdata->eos->getTemperature(pad->pressure, 1./pad->volume);
    }
    cout<<"[LPSolver] ComputeTemperature finished!"<<endl;
}

void LPSolver::updateLocalSpacing(){
    pdata *pad;
    size_t li, lpnum = gdata->particle_data->size();

    for(li=0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        pad->localspacing *= pad->volume/pad->volumeT1;
    }
}

void LPSolver::computeCFLCondition(){
    pdata *pad,*pad2;
    size_t li, lpnum = gdata->particle_data->size();

    double mindt = 100;
    double dt;
    double dx;
    double soundspeed;

    for (li = 0; li<lpnum-1; li++){
        pad = &((*gdata->particle_data)[li]);
        pad2 = &((*gdata->particle_data)[li+1]);
        dx = pad2->x - pad->x;
        soundspeed = pad->soundspeed;
        dt = dx/soundspeed;
        if (dt<mindt) mindt = dt;

        cfldt = mindt * cflcoeff;   
    }
}

bool LPSolver::adjustDtByWriteTimeInterval(){
    cout<<"--------[Time Analysis]----------"<<endl;
    cout<<"[LPSolver] Current time = "<<currenttime<<endl;
    cout<<"[LPSolver] Current timestep = "<<timestep<<endl;
    cout<<"[LPSolver] Current cfldt = "<<cfldt<<endl;
    if(currenttime + cfldt >= nextwritetime){
        cfldt = nextwritetime - currenttime;
        currenttime = nextwritetime;
        nextwritetime += writetimeinterval;
        writestep++;
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
    cout<<"-------------------------------"<<endl;
}
