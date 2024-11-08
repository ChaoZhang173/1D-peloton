#ifndef __BOUNDARY_H__
#define __BOUNDARY_H__
#include <unordered_map>
#include <string>
#include <vector>
#include "eos.h"
class Global_Data;
class Boundary {
public:
    //Constructor
    /// Destructor
	virtual ~Boundary() {}
	
	/**
	 * \brief Get a boundary particle based on a fluid particle      
	 * \param [in] x  The x-coordinate of fluid particle
	 * \param [out] xb  The x-coordinate of boundary particle	
	 */
	virtual void UpdateInflowBoundary(Global_Data* gdata, EOS* m_pEOS, double dt, double m_fInitParticleSpacing){return ;}
    virtual void generateBoundaryParticle(Global_Data *gdata, EOS* m_pEOS, double m_fInitParticleSpacing,double dt,double *mass)  {return ;}  


};


class BoundaryFactory {
public:
	/**
     * \brief Defines a function pointer pointing to a function which creates objects in the Boundary family 
	 */
	typedef Boundary* (*GeoCreateFunc)();
	
	/**
     * \brief   Returns reference to a Singleton object of class BoundaryFactory
	 * 
	 * \param   None 
     *
	 * \return  Reference to a static object of class BoundaryFactory 
	 * 
	 *
	 * Example usage: 
	 * \code
	 *          BoundaryFactory& factory = BoundaryFactory::instance();
	 * \endcode
	 *
	 * \note    This function is implemented based on the lazy Singleton design pattern; 
	 *          therefore, only one BoundaryFactory instance is allowed in each program
	 */
	static BoundaryFactory& instance(); 
	
	/**
     * \brief      Registers (links) the boundary name \e name 
	 *		       with the function \e func for creating objects in the Boundary family
	 *			   
	 *             After registration, \e name can be used as an argument in the createBoundary member function
	 *             for creating objects of the linked type in the Boundary family
	 *	           
	 *  
	 * \param [in] name the boundary name 
	 * \param [in] func the function pointer pointing to the function that creates objects of a specific type
	 *             in the Boundary family
	 * 
	 * \return     None  
	 *
	 * \note       Instead of using this function directly, consider using the BoundaryRegistrar class for
	 *		       the purpose of linking a boundary name and a specific class in the Boundary family. 
	 *             The function is kept public in case one wants to use it directly
	 *		  
	 *
	 */
	void registerBoundary(std::string name, GeoCreateFunc func);
	
	/**
     * \brief      This function creates an object of the class linked to the \name  
	 *
	 * \param [in] name the name linked to a specific class in the Boundary family via the 
	 *				    registrerBoundary member function
	 *  
	 * \return     A Boundary * pointer pointing to an object of a specific class in the Boundary family 
	 * 
	 * Example usage: 
	 * \code
	 *            BoundaryFactory& factory = BoundaryFactory::instance();
	 *            Boundary* newBoundary = factory.createBoundary(name);
	 * \endcode
	 *
	 */
	Boundary* createBoundary(std::string name);
private:
	std::unordered_map<std::string,GeoCreateFunc> bTable; ///< hash table for the (name,creatFunction) pair
	BoundaryFactory() {} ///< for singleton design pattern
	BoundaryFactory(const BoundaryFactory& other); ///< for singleton design pattern (Don't implement)
	BoundaryFactory& operator=(const BoundaryFactory& other); ///< for singleton design pattern (Don't implement)
};


#endif
