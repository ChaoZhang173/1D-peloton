#ifndef __BOUNDARY_PELLET_H__
#define __BOUNDARY_PELLET_H__

#include "boundary.h"
#include "eos.h"
#include <vector>

class PelletInflowBoundary: public Boundary {
public:
    PelletInflowBoundary();	
    virtual ~PelletInflowBoundary() {};
//	virtual void UpdateInflowBoundary(Global_Data* gdata, EOS* m_pEOS, double dt, double m_fInitParticleSpacing);
    virtual void generateBoundaryParticle(Global_Data *gdata, EOS* m_pEOS, double m_fInitParticleSpacing, double dt); 

private:
    double Pinflow;//inflow pressure, constant
	double Uinflow;//inflow velocity, calculated using energy absorb rate
	double Vinflow;//inflow specific volume, constant
};

class PelletOutflowBoundary: public Boundary {

    public:
        PelletOutflowBoundary();
        virtual ~PelletOutflowBoundary() {};
        virtual void generateBoundaryParticle(Global_Data *gdata, EOS* m_pEOS, double m_fInitParticleSpacing, double dt); 
    
    private:
        double thickness;
};

#endif
