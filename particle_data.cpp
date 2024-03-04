#include <iostream>
#include <time.h>
#include <math.h>
#include <algorithm>
#include "particle_data.h"

using namespace std;

Global_Data::Global_Data(Initializer *init) {


}

void Global_Data::initFluidParticles_line(){

  
}

void Global_Data::updateParticleStates(){
    pdata_t *pad;
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