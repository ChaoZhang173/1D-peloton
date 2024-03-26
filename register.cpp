#include "register.h"
#include "state_pellet.h"
#include "state_test.h"
#include "boundary_pellet.h"
namespace{

    StateRegistrar<PelletState> s1("pelletstate");
    StateRegistrar<TestState> s2("teststate");
    BoundaryRegistrar<PelletInflowBoundary> b1("pelletinflowboundary");
}
