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
    // to be finished
    setEOS();

    string statename;
    statename = init->getStateName();
    state= StateFactory::instance().createState(statename);

    // boundary objects
    boundarynumber = init->getBoundaryNumber();
    for(int i = 0; i < boundarynumber; i++){
        string boundaryname = init->getBoundaryName(i);
        boundary.push_back(BoundaryFactory::instance().createBoundary(boundaryname));
    }

}

void Global_Data::initFluidParticles_line(){
    // the total number of particles to be initialized
    int pnum = static_cast<int>(initiallayerlength/initialspacing)+1;
    // Allocate the vector on the heap and assign it to the unique_ptr
    unique_ptr<vector<pdata>> particle_data = make_unique<vector<pdata>>(pnum);
    
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