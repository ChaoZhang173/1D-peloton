#include <iostream>
#include <fstream>
#include "particle_viewer.h"

using namespace std;

ParticleViewer::ParticleViewer(Global_Data *g, const std::string &filename){
    gdata = g;
    outputfilename = filename;
    writestep = 0;
}

void ParticleViewer::writeResult(int step, double time){

    writestep = step;
    string filenameWithStep = outputfilename + "_" + to_string(writestep);
    
    ofstream outputFile(filenameWithStep, ios::app);
    if(!outputFile.is_open()){
        cerr << "[ParticleViewer] Cannot open file " << filenameWithStep << endl;
        return;
    }

    outputFile << "Time: " << time << endl;
    outputFile <<"x\tv\tpressure\tsoundspeed\ttemperature\tvolume\tlocalspacing\n";

    const vector<pdata> &partiles = *(gdata->particle_data);
    for(const pdata &p: partiles){
        outputFile << p.x << "\t" << p.v << "\t" << p.pressure << "\t" << p.soundspeed << "\t" << p.temperature << "\t" << p.volume << "\t" << p.localspacing << endl;
    }

    outputFile.close();
    cout<<"[ParticleViewer] Write file: "<<filenameWithStep<<endl;
} 