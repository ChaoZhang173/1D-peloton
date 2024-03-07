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
    ~Initializer();

    double getCFLcoeff(){return i_CFLCoeff;};
    double getTstart(){return i_StartTime;};
    double getTend(){return i_EndTime;};
    double getWriteTimeInterval(){return i_WriteStep;};

    int getHeatingSourceNumber(){return i_HeatingSourceNumber;};
    double getTeinf(int i){return i_Teinf[i];};
    double getNeinf(int i){return i_Neinf[i];};
    double getWarmupTime(){return i_WarmUpTime;};
    int getEOSChoice(){return i_EOSChoice;};
    int getMaterialChoice(){return i_MaterialChoice;};
    double getLayerLength(){return i_layerlength;};
    double getGamma(){return i_Gamma;};
    int getMaxParticleNumber(){return i_maxParticleNumber;};
    double getInitialSpacing(){return i_Initialspacing;};


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
    void setWriteStep(double w){i_WriteStep = w;};
    void setHeatingSourceNumber(int n){i_HeatingSourceNumber = n;};
    void setWarmupTime(double t){i_WarmUpTime = t;};
    void setEOSChoice(int c){i_EOSChoice = c;};
    void setMaterialChoice(int c){i_MaterialChoice = c;};
    void setTeinf(double t){i_Teinf.push_back(t);};
    void setNeinf(double n){i_Neinf.push_back(n);};
    void setLayerLength(double l){i_layerlength = l;};
    void setGamma(double g){i_Gamma = g;};
    void setMaxParticleNumber(int n){i_maxParticleNumber = n;};
    void setInitialSpacing(double s){i_Initialspacing = s;};

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

    //! a length that used to generate initial particles
    double i_layerlength;
    //! the estimated max number of particles
    int i_maxParticleNumber;
};
#endif
