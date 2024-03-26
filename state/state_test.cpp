#include <iostream>

#include "state_test.h"

TestState::TestState() :
    m_fDen(1. / 197.681), m_fVel(14.5343), m_fPressure(9.37944), xcen(0) {}

double TestState::pressure() {
    return m_fPressure;
}

double TestState::density() {
    return m_fDen;
}

double TestState::velocity() {
    return m_fVel;
}
