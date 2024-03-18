#ifndef __REGISTRAR_H__
#define __REGISTRAR_H__

#include "state.h"
#include "initializer.h"
#include "boundary.h"

template<typename Derived>
class StateRegistrar {
public:
	/**
     * \brief  constructor 
	 *
	 * Registers the \e name given in the argument list.
	 * Registration means linking \e name to a class in the State family.
	 *
	 * \param  name a state name           
	 *
	 * Example usage: Registers \e a_name with class \e SomeState
	 * \code
	 *         StateRegistrar<SomeState> s("a_name");
	 * \endcode
	 */
	StateRegistrar(std::string name);
	/**
     * \brief   This function creates an object of the \e Derived class and returns a pointer to the State class 
	 *
	 * \param   None   
	 *  
	 * \return  A State * pointer that points to an object of \e Derived class  
	 *
	 */
	static State* createFunc();	
};


////////////////////////////////////////////////////////////////////////////////////////
// Start of StateRegistrar
////////////////////////////////////////////////////////////////////////////////////////

template<typename Derived>
State* StateRegistrar<Derived>::createFunc() {
	return new Derived();
}

template<typename Derived>
StateRegistrar<Derived>::StateRegistrar(std::string name) {
	StateFactory& factory = StateFactory::instance(); // note here is the difference from the GeometryRegistrar class
	factory.registerState(name, StateRegistrar<Derived>::createFunc);
}

////////////////////////////////////////////////////////////////////////////////////////
// End of StateRegistrar
////////////////////////////////////////////////////////////////////////////////////////



template<typename Derived>
class BoundaryRegistrar {
public:
	/**
     * \brief  constructor 
	 *
	 * Registers the \e name given in the argument list.
	 * Registration means linking \e name to a class in the Boundary family.
	 *
	 * \param  name a boundary name           
	 *
	 * Example usage: Registers \e a_name with class \e SomeBoundary
	 * \code
	 *         StateRegistrar<SomeBoundary> s("a_name");
	 * \endcode
	 */
	BoundaryRegistrar(std::string name);
	/**
     * \brief   This function creates an object of the \e Derived class and returns a pointer to the Boundary class 
	 *
	 * \param   None   
	 *  
	 * \return  A Boundary * pointer that points to an object of \e Derived class  
	 *
	 */
	 static Boundary* createFunc();

    

};


////////////////////////////////////////////////////////////////////////////////////////
// Start of StateRegistrar
////////////////////////////////////////////////////////////////////////////////////////

template<typename Derived>
Boundary* BoundaryRegistrar<Derived>::createFunc() {
	return new Derived();
}

template<typename Derived>
BoundaryRegistrar<Derived>::BoundaryRegistrar(std::string name) {
	BoundaryFactory& factory = BoundaryFactory::instance(); // note here is the difference from the GeometryRegistrar class
	factory.registerBoundary(name, BoundaryRegistrar<Derived>::createFunc);
}
#endif // __REGISTRAR_H__
