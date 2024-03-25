#include <iostream>
#include <algorithm>
#include <cassert>

#include "particle_data.h"
#include "pellet_solver.h"

using namespace std;

Global_Data::Global_Data(Initializer *init) {

    initiallayerlength = init->getLayerLength();
    initialspacing = init->getInitialSpacing();
    mindx = init->getMinDx();
    //maxparticlenum = init->getMaxParticleNumber();
    gamma = init->getGamma();
    // get background pressure
    backgroundpressure = init->getBackgroundPressure();
    
    // get invalid pressure and density
    invalidpressure = init->getInvalidPressure();
    invaliddensity = init->getInvalidDensity();

    // get the bad pressure and volume, used to delete particles
    badpressure = init->getBadPressure();
    badvolume = init->getBadVolume();

    eoschoice = init->getEOSChoice();
    // set up eos
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
    int pnum = static_cast<int>(initiallayerlength/initialspacing);
    // Allocate the vector on the heap and assign it to the unique_ptr
    particle_data = make_unique<vector<pdata>>(pnum);
    // allocate ghost particle data 
    ghostparticle_data = make_unique<vector<pdata>>(2);
    
    pdata *pad;
    // the first particle is set near but not at the pellet surface
    for (int i = 0; i < pnum; i++){
        pad = &((*particle_data)[i]);
        pad->x = (i+1) * initialspacing;
        pad->v = state->velocity();
        pad->volume = 1./state->density();
        pad->pressure = state->pressure();
        pad->localspacing = initialspacing;
        pad->mass = state->density()*initialspacing;
        pad->soundspeed = eos->getSoundSpeed(pad->pressure, 1./pad->volume, i);
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

    // sort particles based on location
    sort(particle_data->begin(), particle_data->end(), [](const pdata &a, const pdata &b){
        return a.x < b.x;
    });

    // set neighbours
    // left boundary
    pad = &((*particle_data)[0]);
    pad->leftneighbour = &((*ghostparticle_data)[0]);
    pad->rightneighbour = &((*particle_data)[1]);
    // right boundary
    pad = &((*particle_data)[lpnum-1]);
    pad->leftneighbour = &((*particle_data)[lpnum-2]);
    pad->rightneighbour = &((*ghostparticle_data)[1]);
    // other particles
    for(li = 1; li<lpnum-1; li++){
        pad = &((*particle_data)[li]);
        
        pad->leftneighbour = &((*particle_data)[li-1]);
        pad->rightneighbour = &((*particle_data)[li+1]);
    }
}

void Global_Data::generateGhostParticles(){
    pdata *pad,*pad2;
    pdata *ghostpad;
    pellet_info *pellet;

    double dis = 0;

    double vacumm_volume = 1.0e6;

    // get the pellet information, currently only has 1 pellet
    pellet = &((*pellet_solver->pelletlist)[0]);
    
    // left ghost particle
    pad = &((*particle_data)[0]);
    pad2 = &((*particle_data)[1]);
    dis = pad2->x - pad->x;
    ghostpad = &((*ghostparticle_data)[0]);
    ghostpad->x = pad->x - dis;
    ghostpad->localspacing = pad->localspacing;
    ghostpad->ifboundary = true;
    ghostpad->mass = pad->mass;
//    // if the 1st particle is near the pellet surface
//    if(pad->x<10*initialspacing){
        ghostpad->v = pellet->pelletvelocity;
        ghostpad->volume = pellet->vinflow;
        ghostpad->pressure = pellet->pinflow;
        ghostpad->soundspeed = eos->getSoundSpeed(ghostpad->pressure, 1./ghostpad->volume, 0);
//    }
//    // if the 1st particle is far from the pellet surface
//    else{
//        ghostpad->v = pad->v;
//        ghostpad->volume = vacumm_volume;
//        ghostpad->pressure = pad->pressure;
//        ghostpad->soundspeed = pad->soundspeed;
//    }

    // right ghost particle, vacumm
    pad = &((*particle_data)[particle_data->size()-1]);
    pad2 = &((*particle_data)[particle_data->size()-2]);
    dis = pad->x - pad2->x;
    ghostpad = &((*ghostparticle_data)[1]);
    ghostpad->x = pad->x + dis;
    ghostpad->v = pad->v;
    ghostpad->volume = vacumm_volume;
    ghostpad->pressure = backgroundpressure;
    ghostpad->localspacing = pad->localspacing;
    ghostpad->mass = 0.0;
    ghostpad->soundspeed = pad->soundspeed;
    ghostpad->ifboundary = true;


}

void Global_Data::deleteBadParticles(){
    pdata *pad;
    int li;
    size_t lpnum = particle_data->size();

    for(li = lpnum-1; li >= 0; li--){
        pad = &((*particle_data)[li]);
        if(pad->volume > badvolume){
            cout<<"----------Delete Bad Particles----------"<<endl;
            cout<<"Warning: delete particle "<<li<<"at "<<pad->x<<endl;
            cout<<"P = "<<pad->pressure<<"V = "<<pad->volume<<"T = "<<pad->temperature<<endl;
            // delete the particle
            particle_data->erase(particle_data->begin()+li);
            cout<<"-------Finish Delete Bad Particles-------"<<endl;
        }

    }
}

void Global_Data::updatelocalSpacing(){
    pdata *pad;
    size_t li, lpnum = particle_data->size();
    pellet_info *pellet;

    pellet = &((*pellet_solver->pelletlist)[0]);

    // 1st particle
    pad = &((*particle_data)[0]);
    pad->localspacing = 0.5 * ((*particle_data)[1].x - pellet->x);
    // last particle
    pad = &((*particle_data)[lpnum-1]);
    pad->localspacing = pad->x - (*particle_data)[lpnum-2].x;

    for(li=1; li<lpnum-1; li++){
        pad = &((*particle_data)[li]);
        if(pad->x - (*particle_data)[li-1].x < mindx){
            cout<<"[Local Spacing] Warning: re locate particle"<<endl;
            cout<<"[Local Spacing] pad->x = "<<pad->x<<" for particle: "<<li<<endl;
            pad->x = (*particle_data)[li-1].x + mindx;
        } 
        if((*particle_data)[li+1].x - pad->x < mindx){
            cout<<"[Local Spacing] Warning: re locate particle"<<endl;
            cout<<"[Local Spacing] pad->x = "<<pad->x<<" for particle: "<<li<<endl;
            pad->x = (*particle_data)[li+1].x - mindx;
        }
        pad->localspacing = 0.5*((*particle_data)[li+1].x - (*particle_data)[li-1].x);
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

Global_Data::~Global_Data(){
    delete pellet_solver;
    delete eos;
    delete state;
    delete initializer;
    delete boundary[0];
}