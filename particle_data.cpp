#include <iostream>
#include <time.h>
#include <math.h>
#include <algorithm>
#include "particle_data.h"

using namespace std;

Global_Data::Global_Data(Initializer *init) {

    initiallayerlength = init->getLayerLength();
    initialspacing = init->getInitialSpacing();
    //maxparticlenum = init->getMaxParticleNumber();
    gamma = init->getGamma();
    // get background pressure
    backgroundpressure = init->getBackgroundPressure();
    // to be finished
    setEOS();

    // at the first step
    ifStart = true;

    string statename;
    statename = init->getStateName();
    state= StateFactory::instance().createState(statename);

    // boundary objects
    boundarynumber = init->getBoundaryNumber();
    string boundaryname;
    for(int i = 0; i < boundarynumber; i++){
        boundaryname = init->getBoundaryName(i);
        boundary.push_back(BoundaryFactory::instance().createBoundary(boundaryname));
    }

}

void Global_Data::initFluidParticles_line(){
    // the total number of particles to be initialized
    int pnum = static_cast<int>(initiallayerlength/initialspacing)+1;
    // Allocate the vector on the heap and assign it to the unique_ptr
    unique_ptr<vector<pdata>> particle_data = make_unique<vector<pdata>>(pnum);
    // allocate ghost particle data 
    unique_ptr<vector<pdata>> ghostparticle_data = make_unique<vector<pdata>>(2);
    
    pdata *pad;
    // the first particle is set near but not at the pellet surface
    for (int i = 0; i < pnum; i++){
        pad = &((*particle_data)[i]);
        pad->x = i * initialspacing;
        pad->v = state->velocity();
        pad->volume = 1./state->density();
        pad->pressure = state->pressure();
        pad->localspacing = initialspacing;
        pad->mass = state->density()*initialspacing;
        pad->soundspeed = eos->getSoundSpeed(pad->pressure, 1./pad->volume);
        pad->ifboundary = false;
    }
}

void Global_Data::updateParticleStates(){
    pdata *pad;
    int li, lpnum = particle_data->size();

    for(li = 0; li < lpnum; li++){
        pad = &((*particle_data)[li]);
        swap(pad->pressure, pad->pressureT1);
        swap(pad->volume, pad->volumeT1);
        swap(pad->soundspeed, pad->soundspeedT1);
        swap(pad->v, pad->oldv);
    }
}

void Global_Data::reorderParticles(){
    sort(particle_data->begin(), particle_data->end(), [](const pdata &a, const pdata &b){
        return a.x < b.x;
    });
}

void Global_Data::generateGhostParticles(){
    pdata *pad;
    pdata *ghostpad;
    
    double vacumm_volume = 1.0e6;

    // left ghost particle
    pad = &((*particle_data)[0]);
    ghostpad = &((*ghostparticle_data)[0]);
    ghostpad->x = pad->x - 0.5*pad->localspacing;
    ghostpad->v = pad->v;
    ghostpad->volume = pad->volume;
    ghostpad->pressure = pad->pressure;
    ghostpad->localspacing = pad->localspacing;
    ghostpad->mass = pad->mass;
    ghostpad->soundspeed = pad->soundspeed;
    ghostpad->ifboundary = true;

    // right ghost particle, vacumm
    pad = &((*particle_data)[particle_data->size()-1]);
    ghostpad = &((*ghostparticle_data)[1]);
    ghostpad->x = pad->x + 0.5*pad->localspacing;
    ghostpad->v = pad->v;
    ghostpad->volume = vacumm_volume;
    ghostpad->pressure = backgroundpressure;
    ghostpad->localspacing = pad->localspacing;
    ghostpad->mass = 0.0;
    ghostpad->soundspeed = pad->soundspeed;
    ghostpad->ifboundary = true;
    
}