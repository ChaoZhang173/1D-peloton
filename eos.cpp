#include <iostream>
#include <cmath>
#include <cassert>

#include "eos.h"

using namespace std;

double PolytropicEOS::getEnergy(double pressure, double density){
    if(((gamma - 1.)*density)!=0)
        return pressure/((gamma - 1.)*density);
    else {
        cout <<"[EOS] Error: density is zero in PolytropicEOS::getEnergy" << endl;
        cout<<"[EOS] gamma = "<<gamma<<", pressure = "<<pressure<<", density = "<<density<<endl;
        assert(false);
    }
}

double PolytropicEOS::getTemperature(double pressure, double density){
    double R,mu;
    // neon molar mass: g/mol
    mu = 20.18;
    // universal gas constant: in cgs
    R = 83.14;
    // temperature: in eV (from kelvin to ev: 1eV = 11604.525K)
    return mu*pressure/(R*density)/11604.525;
}

// -1: left ghost partilc, -2: right ghost particle
double PolytropicEOS::getSoundSpeed(double pressure, double density, int num){
    double soundspeed;

    if(density != 0)
        soundspeed = gamma*pressure/density;
    else {
        cout<<"[EOS] Error: density is zero in PolytropicEOS::getSoundSpeed"<<endl;
        cout<<"[EOS] gamma = "<<gamma<<", pressure = "<<pressure<<", density = "<<density<<endl;
        cout<<"[EOS] particle number = "<<num<<endl;
        assert(false);
    }
    if(soundspeed > 0)
        return sqrt(soundspeed);
    else if (soundspeed > -1e-10 && soundspeed < 1e-10)
        return 0;
    else {
        cout<<"[EOS] Error: soundspeed is negative in PolytropicEOS::getSoundSpeed"<<endl;
        cout<<"[EOS] gamma = "<<gamma<<", pressure = "<<pressure<<", density = "<<density<<endl;
        cout<<"[EOS] particle number = "<<num<<endl;
        assert(false);
    }
}

double PolytropicEOS::getElectricConductivity(double pressure, double density){
    return 0;
}