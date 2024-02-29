#ifndef __PARTICLE_VIEWER_H__
#define __PARTICLE_VIEWER_H__

#include "particle_data.h"

/*
This file is used to output data.
features:
1. output data file
2. print plot using matlab/pyton
*/

class ParticleViewer {
public:
ParticleViewer(Global_Data * data,const std::string& outputfileName="");

~ParticleViewer(){}

void writeResult(int step,double time);

};
#endif