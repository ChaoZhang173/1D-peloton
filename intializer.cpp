#include <iostream>

#include "initializer.h"

using namespace std;

Initializer::Initializer() {
  setInputs();
  setParams();
}

void Initializer::setParams() {
  
}

void Initializer::readInputfile(const std::string& inputfileName) {
  // to be finished
}

void Initializer::setInputs() {
  setStartTime(0.0);
  setEndTime(0.1);
  setCFLCoeff(0.5);
  setWriteStep(0.05);
  setGamma(1.667);
  setHeatingSourceNumber(1);
  setTeinf(2000);
  setNeinf(1e14);
  setWarmupTime(0.01);
  setEOSChoice(1);
  setMaterialChoice(1);
  setLayerLength(0.01);
  setMaxParticleNumber(10000);
  setInitialSpacing(0.0001);
}