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


    }


}

