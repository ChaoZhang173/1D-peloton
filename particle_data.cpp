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

}

void Global_Data::initFluidParticles_line(){
    // the total number of particles to be initialized
    int pnum = static_cast<int>(initiallayerlength/initialspacing)+1;
    // Allocate the vector on the heap and assign it to the unique_ptr
    unique_ptr<vector<pdata>> particle_data = make_unique<vector<pdata>>(pnum);
    
    pdata *pad;
    // the first particle is set as the pellet surface
    pad = &((*particle_data)[0]);
    pad->x = 0.0;
    pad->v = state->velocity();
    pad->volume = 1./state->density();
    pad->pressure = state->pressure();
    pad->localspacing = initialspacing;
    pad->soundspeed = eos->getSoundSpeed(pad->pressure, 1./pad->volume);
    pad->ifboundary = true;
    
    for (int i = 1; i < pnum; i++){
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
        if (pad->ifboundary) continue;
        swap(pad->pressure, pad->pressureT1);
        swap(pad->volume, pad->volumeT1);
        swap(pad->soundspeed, pad->soundspeedT1);
        swap(pad->v, pad->oldv);
    }
}