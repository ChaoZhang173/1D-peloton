#include <iostream>

#include "state_pellet.h"

PelletState::PelletState():
m_fDen(1./100), m_fVel(0), m_fPressure(16), xcen(0){}

double PelletState::pressure() {
	return m_fPressure;
}

double PelletState::density(){
	return m_fDen;
}

double PelletState::velocity(){
	return m_fVel;
}
