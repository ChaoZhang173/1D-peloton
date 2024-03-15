#include <iostream>
#include <algorithm>
#include <cassert>

#include "particle_data.h"

using namespace std;

Global_Data::Global_Data(Initializer *init) {

    initiallayerlength = init->getLayerLength();
    initialspacing = init->getInitialSpacing();
    mindx = init->getMinDx();
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
    size_t li, lpnum = particle_data->size();

    for(li = 0; li < lpnum; li++){
        pad = &((*particle_data)[li]);
        swap(pad->pressure, pad->pressureT1);
        swap(pad->volume, pad->volumeT1);
        swap(pad->soundspeed, pad->soundspeedT1);
        swap(pad->v, pad->oldv);
    }
}

void Global_Data::reorderParticles(){
    
    pdata *pad;
    size_t li, lpnum = particle_data->size();

    sort(particle_data->begin(), particle_data->end(), [](const pdata &a, const pdata &b){
        return a.x < b.x;
    });

    // left boundary
    pad = &((*particle_data)[0]);
    pad->leftneighbour = &((*ghostparticle_data)[0]);
    pad->rightneighbour = &((*particle_data)[1]);
    // right boundary
    pad = &((*particle_data)[lpnum-1]);
    pad->leftneighbour = &((*particle_data)[lpnum-2]);
    pad->rightneighbour = &((*ghostparticle_data)[1]);

    for(li = 1; li<lpnum-1; li++){
        pad = &((*particle_data)[li]);
        
        pad->leftneighbour = &((*particle_data)[li-1]);
        pad->rightneighbour = &((*particle_data)[li+1]);
       
    }
}

void Global_Data::generateGhostParticles(){
    pdata *pad;
    pdata *ghostpad;
    
    double vacumm_volume = 1.0e6;

    // left ghost particle
    pad = &((*particle_data)[0]);
    ghostpad = &((*ghostparticle_data)[0]);
    ghostpad->x = pad->x - pad->localspacing;
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

void Global_Data::updatelocalSpacing(){
    pdata *pad;
    size_t li, lpnum = particle_data->size();

    // 1st particle
    pad = &((*particle_data)[0]);
    pad->localspacing = (*particle_data)[1].x - pad->x;
    // last particle
    pad = &((*particle_data)[lpnum-1]);
    pad->localspacing = pad->x - (*particle_data)[lpnum-2].x;

    for(li=1; li<lpnum-1; li++){
        pad = &((*particle_data)[li]);
        if(pad->x - (*particle_data)[li-1].x < mindx){
            cout<<"[Local Spacing] Warning: replacing particle"<<endl;
            cout<<"[Local Spacing] pad->x = "<<pad->x<<" for particle: "<<li<<endl;
            pad->x = (*particle_data)[li-1].x + mindx;
        } 
        if((*particle_data)[li+1].x - pad->x < mindx){
            cout<<"[Local Spacing] Warning: replacing particle"<<endl;
            cout<<"[Local Spacing] pad->x = "<<pad->x<<" for particle: "<<li<<endl;
            pad->x = (*particle_data)[li+1].x - mindx;
        }
        pad->localspacing = 0.5*((*particle_data)[li-1].x + (*particle_data)[li+1].x);
        if(pad->localspacing < 3*mindx){
            cout<<"[Local Spacing] Warning: localspacing too small!"<<endl;
            cout<<"[Local Spacing] pad->localspacing = "<<pad->localspacing<<endl;
            pad->x = 0.5*((*particle_data)[li-1].x + (*particle_data)[li+1].x);
            cout<<"[Local Spacing] pad->x = "<<pad->x<<" for particle: "<<li<<endl;
        }
    }
}

void Global_Data::setEOS(){
    // 1: polytropic; 2: saha neon;
    if(eoschoice == 1){
        eos = new PolytropicEOS(gamma);
    }
}