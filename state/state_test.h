#ifndef __STATE_TEST_H__
#define __STATE_TEST_H__

#include "state.h"
#include <math.h>

/* 
This file is for initializing the test state.
*/

class TestState: public State{
public:
	TestState();
	virtual ~TestState(){};
	virtual double pressure();
	virtual double density();
	virtual double velocity();
private:
	double m_fDen;
	double m_fVel;
	double m_fPressure;
	double xcen;
};
#endif

