#include <iostream>

#include "initializer.h"

using namespace std;

Initializer::Initializer() {
  setInputs();
  setParams();
}

void Initializer::setParams() {
  return ;
}

void Initializer::readInputfile(const std::string& inputfileName) {
  // to be finished
}

void Initializer::setInputs() {
  setStartTime(0.0);
  setEndTime(0.1);
  setCFLCoeff(0.5);
  // set the time interval to write data
  setWriteStep(0.05);
  setGamma(1.667);
  // set the number of heating sources
  setHeatingSourceNumber(1);
  // set the temperature and density of heating sources
  setTeinf(2000);
  setNeinf(1e14);
  // set the time to warm up the system, used in the heatingmodel
  setWarmupTime(0.01);
  // set the choice of EOS, 1 for polytropic, 2 for saha neon,
  setEOSChoice(1);
  // set pellet material; 0:neon
  setMaterialChoice(0);
  // the length of layer that will generate initial particles
  setLayerLength(0.01);
  // set the estimated max number of particles, currenlyt not used
  //setMaxParticleNumber(10000);
  // set the initial spacing between particles
  setInitialSpacing(0.0001);
  // set the smallest dx between particles
  setMinDx(1e-6);
  setBackgroundPressure(0.64);//2kev 1e14 -> 0.64
  // set state name
  setStateName("pelletstate");
  setBoundaryNumber(1);// currently only one boundary
  setBoundaryName("pelletinflowboundary");
  // set pellets, currently only use 1 pellet
  setPellentNumber(1);
}