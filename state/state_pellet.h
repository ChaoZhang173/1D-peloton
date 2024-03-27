#ifndef __STATE_PELLET_H__
#define __STATE_PELLET_H__

#include "state.h"
#include <math.h>

/* 
This file is for initializing the pellet state.
*/

class PelletState: public State{
public:
	PelletState();
	virtual ~PelletState(){};
	virtual double pressure();
	virtual double density();
	virtual double velocity();
private:
	double m_fDen;
	double m_fVel;
	double m_fPressure;
	double m_xcen;
};
#endif

