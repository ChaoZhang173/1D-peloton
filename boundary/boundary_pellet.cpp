#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "boundary_pellet.h"
#include "particle_data.h"
#include "pellet_solver.h"

using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(16),Uinflow(0),Vinflow(100){}

void PelletInflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx, double dt){

    PelletSolver *p = g->pellet_solver;
    pdata *pad;
    pellet_info *pellet;

    int pi, pnum = p->pelletlist->size();
    int li, lpnum = g->particle_data->size();

    double mass_fix = dx/Vinflow;
    double massflowrate;
    double pelletvinflow, pelletpinflow;

    // compute massflowrate, need to be finished
    p->computeMassFlowRate();
    // currently only has 1 pellet
    for(pi = 0; pi < pnum; pi++){
        pellet = &((*p->pelletlist)[pi]);
        massflowrate = pellet->massflowrate;
        // if not at the beggining of the job, and massflowrate >0, pinflow is a number  
        if(!g->ifStart&& abs(massflowrate)>1e-10 &&!isnan(pellet->pinflow)){
            pelletvinflow = pellet->vinflow;
            pelletpinflow = pellet->pinflow;
        } 
        else{ // if at the beggining of the job, or massflowrate <=0, pinflow is not a number
            pellet->pelletvelocity = 0.;
            pellet->vinflow = pelletvinflow = Vinflow;
            pellet->pinflow = pelletpinflow = Pinflow;
        }
        
        // output the particle_data[0]->x
        cout<<"[Boundary] the smallest x = "<<g->particle_data->at(0).x<<endl;
        // generate new particles
        int newParticleNum = (int)(massflowrate*dt/mass_fix);
        double dis = pellet->pelletvelocity*dt; // layer that will generate new particles
        cout<<"[Boundary] layer thickness = "<<dis<<endl;
        cout<<"[Boundary] new paticle number = "<<newParticleNum<<endl;

        double actualdx = max(pelletvinflow*mass_fix,0.0);

        // access the vector inside the unique_ptr
        vector<pdata> *particle_data = g->particle_data.get();
        // currently generate particle at 0~smallest x
        double newdx = g->particle_data->at(0).x/newParticleNum;
        for(li=0; li<newParticleNum; li++){
            pad->x = li*newdx;
            pad->v = pellet->pelletvelocity;
            pad->volume = pelletvinflow;
            pad->pressure = pelletpinflow;
            pad->localspacing = actualdx;
            pad->mass = mass_fix;
            pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure, 1./pad->volume);
            pad->ifboundary = false;
            particle_data->emplace_back(*pad);
        }
    }

}

