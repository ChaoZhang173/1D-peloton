#ifndef __INITIALIZER_H__
#define __INITIALIZER_H__

#include <vector>

using namespace std;
class Global_Data;

/*
This file is for initializer class. 
It is used to initialize the parameters and pass them to Global_Data class.
In later version may use a input file to initialize, but currently only set to default values.

The variables in private member are starts with i_ to indicate that they are initialized by initializer.
*/

class Initializer {
    friend class Global_Data; // Global_Data has access to private,public and proteced members
    public:
    // constructor
    Initializer();
    // destructor
    ~Initializer() {};

    double getCFLcoeff(){return i_CFLCoeff;};
    double getTstart(){return i_StartTime;};
    double getTend(){return i_EndTime;};
    double getWriteTimeInterval(){return i_WriteStep;};

    int getHeatingSourceNumber(){return i_HeatingSourceNumber;};
    double getTeinf(int i){return i_Teinf[i];};
    double getNeinf(int i){return i_Neinf[i];};
    double getWarmupTime(){return i_WarmUpTime;};
    int getEOSChoice(){return i_EOSChoice;};
    //! pellet material; 0:neon 
    int getMaterialChoice(){return i_MaterialChoice;};
    double getLayerLength(){return i_layerlength;};
    double getGamma(){return i_Gamma;};
    int getMaxParticleNumber(){return i_maxParticleNumber;};
    double getInitialSpacing(){return i_Initialspacing;};
    string getStateName(){return i_StateName;};
    int getBoundaryNumber(){return i_BoundaryNumber;};
    string getBoundaryName(int i) {return i_BoundaryName[i];}

    int getPelletNumber(){return i_PelletNumber;};
    double getPelletLocation(){return i_PelletLocation;};
    
    double getBackgroundPressure(){return i_BackgroundPressure;};
    double getMinDx(){return i_minDx;};

    double getInvalidPressure(){return i_invalidPressure;};
    double getInvalidDensity(){return i_invalidDensity;};

    double getBadPressure(){return i_badPressure;};
    double getBadVolume(){return i_badVolume;};

    private:
    // currently not using
    void readInputfile(const std::string& inputfileName);
    //! function that used to set some initial values, take replace of input file
    void setInputs();
    // ! function that used to set some deafult values that are not from input
    void setParams();

    // functions that used to set vaules (take replace of input file)
    void setStartTime(double t){i_StartTime = t;};
    void setEndTime(double t){i_EndTime = t;};
    void setCFLCoeff(double c){i_CFLCoeff = c;};
    //! set the time interval to write data
    void setWriteStep(double w){i_WriteStep = w;};
    //! how many heating sources
    void setHeatingSourceNumber(int n){i_HeatingSourceNumber = n;};
    //! the time to warm up the system, used in the heatingmodel
    void setWarmupTime(double t){i_WarmUpTime = t;};
    //! the choice of EOS, 1 for polytropic, 2 for saha neon,
    // currently only use 1
    void setEOSChoice(int c){i_EOSChoice = c;};
    //! pellet material; 0:neon 
    void setMaterialChoice(int c){i_MaterialChoice = c;};
    //! heating source temperature (electrons)
    void setTeinf(double t){i_Teinf.push_back(t);};
    //! heating source density (electrons)
    void setNeinf(double n){i_Neinf.push_back(n);};
    //! the length of layer that will generate initial particles
    void setLayerLength(double l){i_layerlength = l;};
    void setGamma(double g){i_Gamma = g;};
    //! the estimated max number of particles, currently not used
    void setMaxParticleNumber(int n){i_maxParticleNumber = n;};
    //! the initial distance between particles
    void setInitialSpacing(double s){i_Initialspacing = s;};
    //! the state name, currently only use "pelletstate"
    void setStateName(string s){i_StateName = s;};
    //! the number of boundary objects, currently only use 1
    void setBoundaryNumber(int n){i_BoundaryNumber = n;};
    //! the name of boundary object, currently only use "pelletinflowboundary"
    void setBoundaryName(string s){i_BoundaryName.push_back(s);};
    //! the number of pellets, currently only use 1
    void setPelletNumber(int n){i_PelletNumber = n;};
    //! the location of pellet, currently only use 0
    void setPelletLocation(double l){i_PelletLocation = l;};
    //! the background pressure
    void setBackgroundPressure(double p){i_BackgroundPressure = p;};
    //! the smallest dx between particles
    void setMinDx(double dx){i_minDx = dx;};
    //! set the invalid pressure
    void setInvalidPressure(double p){i_invalidPressure = p;};
    //! set the invalid density
    void setInvalidDensity(double d){i_invalidDensity = d;};
    //! set the bad pressure
    void setBadPressure(double p){i_badPressure = p;};
    //! set the bad volume
    void setBadVolume(double v){i_badVolume = v;};

    double i_StartTime;
    double i_EndTime;
    double i_CFLCoeff;
    double i_WriteStep;//! the time interval to write data
    double i_Gamma;
    double i_Initialspacing;
    int i_HeatingSourceNumber;
    double i_WarmUpTime;
    int i_EOSChoice;
    int i_MaterialChoice;
    vector<double> i_Teinf;
    vector<double> i_Neinf;
    double i_BackgroundPressure;

    double i_invalidPressure;
    double i_invalidDensity;

    double i_badPressure;
    double i_badVolume;

    double i_minDx;

    string i_StateName;

    int i_BoundaryNumber;
    vector<string> i_BoundaryName;

    //! number of pellets
    int i_PelletNumber;
    //! the location of pellet
    double i_PelletLocation;

    //! a length that used to generate initial particles
    double i_layerlength;
    //! the estimated max number of particles
    int i_maxParticleNumber;
};
#endif
