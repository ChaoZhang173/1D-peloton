#ifndef __EOS_H__
#define __EOS_H__

/*
This file is for equation of state class. 
We use virtual function so that they could be overrided by the derived class.
Other classes inherit from the EOS calss and override the virtual functions, they have access
to the protected and public members of the EOS class.
It has the following eos choices:(EOSChoice)
  1. polytropic (currently working)
  2. saha neon (in the future)
  3. others
*/
class EOS {

protected:
  double gamma;
  // the choice of EOS, 1 for polytropic, 2 for saha neon
  int EOSChoice; 
  
public:
  // destructor
  virtual ~EOS(){}; 

  int getEOSChoice() const {return EOSChoice;};
  double getGamma() const {return gamma;};

  virtual double getEnergy(double pressure, double density) = 0;
  virtual double getTemperature(double pressure, double density) = 0;
  virtual double getSoundSpeed(double pressure, double density, int num) = 0;
  virtual double getElectricConductivity(double pressure, double density) = 0;

};

class PolytropicEOS : public EOS {
protected:

public:
  PolytropicEOS(double gamma){EOSChoice = 1; this->gamma = gamma;};

  virtual ~PolytropicEOS(){};
  virtual double getEnergy(double pressure, double density);
  virtual double getTemperature(double pressure, double density);
  virtual double getSoundSpeed(double pressure, double density, int num);
  virtual double getElectricConductivity(double pressure, double density);
};

#endif
