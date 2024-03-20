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
    pdata *pad = new pdata;
    pellet_info *pellet;

    size_t pi;
    size_t pnum = p->pelletlist->size();
    size_t li;

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
        
        cout<<"-----------------"<<endl;
        // output the particle_data[0]->x
        cout<<"[Boundary] the smallest x = "<<g->particle_data->at(0).x<<endl;
        // generate new particles
        size_t newParticleNum = (size_t)(massflowrate*dt/mass_fix);
        cout<<"[Boundary] massflowrate = "<<massflowrate<<endl;
        cout<<"[Boundary] dt = "<<dt<<endl;
        cout<<"[Boundary] mass_fix = "<<mass_fix<<endl;
        cout<<"[Boundary] newParticleNum = "<<newParticleNum<<endl;
        double dis = pellet->pelletvelocity*dt; // layer that will generate new particles
        cout<<"[Boundary] pellet velocity = "<<pellet->pelletvelocity<<endl;
        cout<<"[Boundary] new layer thickness = "<<dis<<endl;
        cout<<"[Boundary] new paticle number = "<<newParticleNum<<endl;
        
        // if no new particles will be generated, return
        if(dis == 0 || newParticleNum == 0){
            cout<<"[Boundary] no particles will be generated"<<endl;
            return;
        }
        cout<<"-----------------"<<endl;
        // used as localspacing in 3D code, not used in 1D code
        //double actualdx = max(pelletvinflow*mass_fix,0.0);

        // access the vector inside the unique_ptr
        vector<pdata> *particle_data = g->particle_data.get();
        // currently generate particle at 1*newdx~pellet velocity *dt
        double newdx = dis/newParticleNum;
        for(li=0; li<newParticleNum; li++){
            pad->x = (li+1)*newdx;
            pad->v = pellet->pelletvelocity;
            pad->volume = pelletvinflow;
            pad->pressure = pelletpinflow;
            pad->localspacing = newdx; // some of them need to be updated
            pad->mass = mass_fix;
            pad->soundspeed = m_pEOS->getSoundSpeed(pad->pressure, 1./pad->volume, li);
            pad->ifboundary = false;
            particle_data->emplace_back(*pad);
        }
    }
    delete pad;
}

