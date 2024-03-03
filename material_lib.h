#ifndef __MATERIAL_LIB__
#define __MATERIAL_LIB__

/*
This file is for material class, currently only have neon, 
will add more materials in the future
*/

class Material{

    public:
        double getMu() const {return mu;} 
        double getMass() const {return mass;}
        double getZ()const {return Z;}
        double getI() const{return I;}
        double getSublimationEnergy() const {return sublimationEnergy;}
        virtual double  getOne_Plus_Zstar(double teinf ) = 0;
    protected: 
        double mu;
        double mass;
        double Z;
        double I;
        double sublimationEnergy;
        double one_plus_Zstar;
};

class Neon: public Material{
    
    public:  
         virtual double  getOne_Plus_Zstar(double teinf ); 
        Neon();
    
        ~Neon();


};



#endif
