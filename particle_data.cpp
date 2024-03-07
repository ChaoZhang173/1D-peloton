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

}

void Global_Data::initFluidParticles_line(){
    // the total number of particles to be initialized
    int pnum = static_cast<int>(initiallayerlength/initialspacing);
    // Allocate the vector on the heap and assign it to the unique_ptr
    unique_ptr<vector<pdata>> particle_data = make_unique<vector<pdata>>();

    pdata pad; 
    // the first particle is set as the pellet surface
    pad.x = 0.0;

  
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