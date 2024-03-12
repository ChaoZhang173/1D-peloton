#include <iostream>
#include <cmath>
#include <cassert>

#include "pellet_solver.h"
#include "geometry_pellet.h"

using namespace std;

PelletSolver::PelletSolver(Initializer *init,Global_Data*g){
    gdata = g;
    setPelletMaterial(init->getMaterialChoice());
    heatsource_numer = init->getHeatingSourceNumber();
    for(int i = 0; i < heatsource_numer; i++){
        teinf.push_back(init->getTeinf(i));
        neinf.push_back(init->getNeinf(i));
    }
    warmuptime = init->getWarmupTime();

    // initialize the pellets
    pelletnumber = init->getPelletNumber();
    pelletlist = make_unique<vector<pellet_info>>(pelletnumber);
    pellet_info *pellet;
    for (int i = 0; i < pelletnumber; i++){
        pellet = &((*pelletlist)[i]);
        pellet->x = -0.5*init->getInitialSpacing();
        pellet->layerlength = init->getLayerLength();
    }
    
}

void PelletSolver::heatingModel(double currenttime){
    cout << "[Heat] Enter heatingModel..." << endl;

    // coumpute density integral
    computeDensityIntegral();
    // compute heat deposition 
    computeHeatDeposition(currenttime);

}

void PelletSolver::computeDensityIntegral(){
    cout << "[Heat] Enter computeDensityIntegral..." << endl;
    
    pdata *pad;
    pdata *pad2; // the right neighbour particle
    int li, lpnum = gdata->particle_data->size();
    double integral_sum = 0.;
    double integral_local = 0.;
    // from right end to the pellet 
    for(li = lpnum-1; li >= 0; li--){
        pad = &((*gdata->particle_data)[li]);
        pad2 = &((*pad->neighbourparticle)[1]);
        integral_local = 0.5*(pad2->x - pad->x)*(1./pad2->volume + 1./pad->volume);
        integral_sum += integral_local;
        pad->leftintegral = integral_sum;

        // the following is set for computeHeatDeposition
        pad->deltaq = 0.;
        pad->qplusminus = 0.;

    }

    cout << "[Heat] Exit computeDensityIntegral..." << endl;
}

void PelletSolver::computeHeatDeposition(double currenttime){
    cout << "[Heat] Enter computeHeatDeposition..." << endl;

    for(int i = 0; i < heatsource_numer; i++){
        cout << "Heating source " << i+1 << " : " << teinf[i] << "ev " << neinf[i] << endl;
        addHeatSource(teinf[i], neinf[i], currenttime);
    }


    cout << "[Heat] Exit computeHeatDeposition..." << endl;
}

void PelletSolver::addHeatSource(double teinf, double neinf,double currenttime){

    pdata *pad;
    int li, lpnum = gdata->particle_data->size();

    // get 1+Z*
    one_plus_Zstar = material->getOne_Plus_Zstar(teinf);
    
    const double e = heatK*(2.99792458e7)/100;
    const double lnLambda = log(2*teinf/I*sqrt(exp(1)/2));
    double tauinf, taueff, tauleft, uleft, qinf, guleft;
    double nt;

    double k_warmup;

    if(currenttime > warmuptime)
        k_warmup = 1.0;
    else 
        k_warmup = currenttime/warmuptime;

    for(li = 0; li < lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        tauleft = pad->leftintegral/mass*Z;
        tauinf = heatK*heatK*teinf*teinf/(8.0*3.1416*e*e*e*e*lnLambda);
        taueff = tauinf/(0.625+0.55*sqrt(one_plus_Zstar));
        uleft = tauleft/taueff;
        qinf = sqrt(2.0/3.1416/masse)*neeff*pow(heatK*teinf,1.5);
        guleft = sqrt(uleft)*Bessel_K1(sqrt(uleft))/4;
        nt = 1.0/pad->volume/mass;

        pad->deltaq += qinf*nt*Z/taueff*guleft*k_warmup;
        pad->qplusminus += qinf*0.5*uleft*Bessel_Kn(2,sqrt(uleft))*k_warmup;
    }
    


}

void PelletSolver::setPelletMaterial(int id){
    materialid = id;
    if(id == 0){
        material = new Neon();
    }
    else{
        assert(false);
    }
    mu = material->getMu();
    mass = material->getMass();
    Z = material->getZ();
    I = material->getI();
    sublimationenergy = material->getSublimationEnergy();
}

void PelletSolver::computeMassFlowRate(){
    
    pellet_info *pellet;
    pdata *pad;
    
    int pi, pnum = pelletlist->size();
    double massflowrate;
    double qsum;
    
    pad = &((*gdata->particle_data)[0]);
    for(pi = 0; pi < pnum; pi++){
        pellet = &((*pelletlist)[pi]);
        if(gdata->ifStart){
            pellet->massflowrate = 0.;
            continue;
        }
        qsum = pad->qplusminus;
        massflowrate = qsum/sublimationenergy;
        pellet->massflowrate = massflowrate;
        cout<<"[Pellet] The mass flow rate = "<<massflowrate<<endl;
    }
}

void PelletSolver::computeBoundaryCondition(Global_Data *g, double dt, double dx){
    
    pdata *pad;
    pellet_info *pellet;
    
    int li, lpnum = g->particle_data->size();
    int pi, pnum = pelletlist->size();

    int counter_nei = 0;
    double qsum = 0;
    double qsum_bc = 0;
    double vol = 0;
    double vol_bc = 0;
    double ur = 0;
    double ur_bc = 0;
    double pres = 0;
    double pres_bc = 0;
    double soundspeed = 0;
    double soundspeed_bc = 0;
    double massflowrate;

    double B,C; 

    // states from charactieristic
    double pvolume, ppressure,pur;

    double v,ss; 
    double sound;
    double r_shift;
    double x; // the position of the particle
    double pellet_cen; 

    double gamma = g->gamma;
    double R = 83.1466/mu;
    double Ts = 450;

    for (pi = 0; pi<pnum; pi++){
        pellet = &((*pelletlist)[pi]);
        counter_nei = 0;
        qsum = 0;
        vol = 0;
        ur = 0;
        pres = 0;
        soundspeed = 0;
        pad = &((*g->particle_data)[0]);
        pellet_cen = -0.5*pad->localspacing;

        for(li = 0; li < lpnum; li++){
            pad = &((*g->particle_data)[li]);
            x = pad->x;
            v = (pad->v+pad->oldv)/2;
            sound = (pad->soundspeed+pad->soundspeedT1)/2;
            if(v == 0)
                ss = 0;
            else{
                ss = sound;
            }
            r_shift = x-pellet_cen - ss*dt;
            if(r_shift<pellet_cen+dx && r_shift>pellet_cen-dx){
                pvolume = pad->volume;
                ppressure = pad->pressure;
                pur = v;
                qsum += computeQplusminuisGradient(pad,pellet);
                vol += pvolume;
                pres += ppressure;
                soundspeed += pad->soundspeed;
                ur += pur;
                counter_nei++;
            }
        }
        if(counter_nei == 0){
            cout<<"[Boundary] No neighbour particles found for pellet "<<pi<<"!"<<endl;
        }
        qsum_bc = qsum/counter_nei;
        vol_bc = vol/counter_nei;
        pres_bc = pres/counter_nei;
        soundspeed_bc = soundspeed/counter_nei;
        ur_bc = ur/counter_nei;
        massflowrate = pellet->massflowrate;
        if(abs(massflowrate) <1e-10){
            B = C = 0;
        }
        else{
            B = (pres_bc+dt*(gamma-1)*(qsum_bc))*vol_bc/soundspeed_bc - ur_bc;
            C = -massflowrate*R*Ts*vol_bc/soundspeed_bc;
        }
        pellet->pelletvelocity = (-B+sqrt(B*B-4*C))/2;
        pellet->vinflow = pellet->pelletvelocity/massflowrate;
        pellet->pinflow = R*Ts/pellet->vinflow;

        cout<<"[Boundary] Inflow volume = "<<vol_bc<<endl;
        cout<<"[Boundary] Inflow pressure = "<<pres_bc<<endl;
        cout<<"[Boundary] Inflow soundspeed = "<<soundspeed_bc<<endl;
        cout<<"[Boundary] Inflow velocity = "<<ur_bc<<endl;
        cout<<"[Boundary] Massflowrate = "<<massflowrate<<endl;
        cout<<"[Boundary] Ts = "<<Ts<<endl;
        cout<<"[Boundary] pellet velcoity = "<<pellet->pelletvelocity<<endl;
    }
}

double PelletSolver::computeQplusminuisGradient(pdata *pad, pellet_info *pellet){
    
}

double Bessel_I0(double x)
{
    double   p1 = 1.0;
    double   p2 = 3.5156229;
    double   p3 = 3.0899424;
    double   p4 = 1.2067492;
    double   p5 = 0.2659732;
    double   p6 = 0.360768e-1;
    double   p7 = 0.45813e-2;

    double   q1 = 0.39894228;
    double   q2 = 0.1328592e-1;
    double   q3 = 0.225319e-2;
    double   q4 = -0.157565e-2;
    double   q5 = 0.916281e-2;
    double   q6 = -0.2057706e-1;
    double   q7 = 0.2635537e-1;
    double   q8 = -0.1647633e-1;
    double   q9 = 0.392377e-2;

    double   ax, y, value;

    if (fabs(x) < 3.75)
    {
        y = (x/3.75)*(x/3.75);//sqr
        value = p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))));
    }
    else
    {
        ax = fabs(x);
        y = 3.75/ax;

        value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
    }

    return value;
}

/* Bessel_I1 returns the modifies Bessel function I1(x) of positive real x  */
double Bessel_I1(double x)
{
    double   p1 = 0.5;
    double   p2 = 0.87890594;
    double   p3 = 0.51498869;
    double   p4 = 0.15084934;
    double   p5 = 0.2658733e-1;
    double   p6 = 0.301532e-2;
    double   p7 = 0.32411e-3;

    double   q1 = 0.39894228;
    double   q2 = -0.3988024e-1;
    double   q3 = -0.362018e-2;
    double   q4 = 0.163801e-2;
    double   q5 = -0.1031555e-1;
    double   q6 = 0.2282967e-1;
    double   q7 = -0.2895312e-1;
    double   q8 = 0.1787654e-1;
    double   q9 = -0.420059e-2;

    double   ax, y, value;

    if (fabs(x) < 3.75)
    {
        y = (x/3.75)*(x/3.75);//sqr
        value = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
    }
    else
    {
        ax = fabs(x);
        y = 3.75/ax;

        value = (exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))));
        if (x < 0)
            value *= -1.0;
    }
    return value;
}

/* Bessel_K0 returns the modifies Bessel function K0(x) of positive real x  */
double Bessel_K0(double x)
{
    double   p1 = -0.57721566;
    double   p2 = 0.4227842;
    double   p3 = 0.23069756;
    double   p4 = 0.348859e-1;
    double   p5 = 0.262698e-2;
    double   p6 = 0.1075e-3;
    double   p7 = 0.74e-5;

    double   q1 = 1.25331414;
    double   q2 = -0.7832358e-1;
    double   q3 = 0.2189568e-1;
    double   q4 = -0.1062446e-1;
    double   q5 = 0.587872e-2;
    double   q6 = -0.25154e-2;
    double   q7 = 0.53208e-3;

    double   y, value;

    if (x <= 2.0)
    {
        y = x*x/4.0;
        value = (-log(x/2.0)*Bessel_I0(x))+(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
    }
    else
    {
        y = 2.0/x;
        value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
    }
    return value;
}

/* Bessel_K1 returns the modifies Bessel function K1(x) of positive real x  */
double Bessel_K1(double x)
{
    double   p1 = 1.0;
    double   p2 = 0.15443144;
    double   p3 = -0.67278579;
    double   p4 = -0.18156897;
    double   p5 = -0.01919402;
    double   p6 = -0.110404e-2;
    double   p7 = -0.4686e-4;

    double   q1 = 1.25331414;
    double   q2 = 0.23498619;
    double   q3 = -0.3655620e-1;
    double   q4 = 0.1504268e-1;
    double   q5 = -0.780353e-2;
    double   q6 = 0.325614e-2;
    double   q7 = -0.68245e-3;

    double   y, value;

    if (x <= 2.0)
    {
        y = x*x/4.0;
        value = (log(x/2.0)*Bessel_I1(x))+(1.0/x)*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))));
    }
    else
    {
        y = 2.0/x;
        value = (exp(-x)/sqrt(x))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*q7))))));
    }
    return value;
}

/* Bessel_Kn returns the modifies Bessel function Kn(x) of positive real x for n >= 2 */
double Bessel_Kn(int n,double x)
{

    int    j;
    double bk, bkm, bkp, tox;

    if (n < 2)
    {
        printf("Error in Bessel_Kn(), the order n < 2: n = %d\n",n);
        assert(false);
        return 0;
    }

    tox = 2.0/x;
    bkm = Bessel_K0(x);
    bk = Bessel_K1(x);

    for (j = 1; j < n; j++)
    {
        bkp = bkm + j*tox*bk;
        bkm = bk;
        bk = bkp;
    }

    return bk;
}