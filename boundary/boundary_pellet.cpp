#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "boundary_pellet.h"
#include "particle_data.h"
#include "pellet_solver.h"

using namespace std;

PelletInflowBoundary::PelletInflowBoundary():Pinflow(16),Uinflow(0),Vinflow(100){}

void PelletInflowBoundary::generateBoundaryParticle(Global_Data *g, EOS* m_pEOS, double dx, double dt, double *mass){

    PelletSolver *p = g->pellet_solver;
    pdata *pad = new pdata;
    pellet_info *pellet;

    double pellet_cen, pellet_radius; 

    size_t pi;
    size_t pnum = p->pelletlist->size();
    size_t li;

    // accumulate mass from last time step
    double *mass_accum = mass;
    double mass_fix = dx/Vinflow;
    double massflowrate;
    double pelletvinflow, pelletpinflow;

    // compute massflowrate, need to be finished
    p->computeMassFlowRate();

    // currently only has 1 pellet
    for(pi = 0; pi < pnum; pi++){
        pellet = &((*p->pelletlist)[pi]);
        pellet_cen = pellet->x;
        pellet_radius = pellet->radius;
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

        // set the boundary condition to fixed value, do test
        pellet->pinflow = 9.37944;
        pellet->vinflow = 197.681;
        pellet->pelletvelocity = 14.5343;
        pellet->massflowrate = 0.03679986557478038;
        massflowrate = pellet->massflowrate;
        pelletpinflow = pellet->pinflow;
        pelletvinflow = pellet->vinflow;
        cout<<"--------[Boundary Test]---------"<<endl;
        cout<<"The pinflow has been set to: "<<pelletpinflow<<endl;
        cout<<"The vinflow has been set to: "<<pelletvinflow<<endl;
        cout<<"The pellet velocity has been set to: "<<pellet->pelletvelocity<<endl;
        cout<<"The massflowrate has been set to: "<<massflowrate<<endl;
        
        cout<<"-----------------"<<endl;
        // output the particle_data[0]->x
        cout<<"[Boundary] the smallest x = "<<g->particle_data->at(0).x<<endl;
        // generate new particles, new mass + accumulated mass
        size_t newParticleNum = (size_t)((massflowrate*dt+(*mass_accum))/mass_fix);
        cout<<"[Boundary] input accumulated mass = "<<*mass_accum<<endl;
        // compute accumulated mass for next time step
        *mass_accum = massflowrate*dt+(*mass_accum) - newParticleNum*mass_fix;
        cout<<"[Boundary] new accumulated mass = "<<*mass_accum<<endl;
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
            pad->x = pellet_cen + pellet_radius + (li+1)*newdx;
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

