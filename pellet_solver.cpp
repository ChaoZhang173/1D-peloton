#include <iostream>
#include <cmath>
#include <cassert>

#include "pellet_solver.h"

using namespace std;

double Bessel_Kn(int n, double x);
double Bessel_K1(double x);
double Bessel_I0(double x);
double Bessel_I1(double x);

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
    size_t li, lpnum = gdata->particle_data->size();
    double integral_sum = 0.;
    double integral_local = 0.;
    // from right end to the pellet 
    for(li = lpnum-1; li >= 0; li--){
        pad = &((*gdata->particle_data)[li]);
        pad2 = pad->rightneighbour;
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
    size_t li, lpnum = gdata->particle_data->size();
    
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
    one_plus_Zstar = material->getOne_Plus_Zstar(teinf[0]);
}

void PelletSolver::computeMassFlowRate(){
    
    pellet_info *pellet;
    pdata *pad;
    
    size_t pi, pnum = pelletlist->size();
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
    
    size_t li, lpnum = g->particle_data->size();
    size_t pi, pnum = pelletlist->size();

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
                qsum += computeQplusminuisGradient(pad);
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

double PelletSolver::computeQplusminuisGradient(pdata *pad){
    
    pdata *pad1, *pad2;
    double dq;

    pad1 = pad->leftneighbour;
    pad2 = pad->rightneighbour;

    dq = (pad2->qplusminus - pad1->qplusminus)/(pad2->x - pad1->x);
    return dq;
}

void PelletSolver::neonRadiationCooling(double dt){
    
    if(gdata->eoschoice != 2){
        cout<<"[Radiation] The neonRadiationCooling is not used!"<<endl;
        return;
    }
    size_t li, lpnum = gdata->particle_data->size();
    double T, rad_cool;
    double sc, pressure, density;
    int timestep = 1;
    pdata *pad;
    for(li = 0; li<lpnum; li++){
        pad = &((*gdata->particle_data)[li]);
        density = 1./pad->volume;
        for(int i=0;i<timestep;i++){
            sc = pad->soundspeed;
            pressure = pad->pressure;
            T = gdata->eos->getTemperature(pressure,density);
            rad_cool = neon_radiation_power_density(density, T);
            pad->radcool = rad_cool;
            pad->pressure = pressure - (sc*sc*density/pressure -1 )*rad_cool*dt/timestep;
            if(pad->pressure < 0){
                pad->pressure = pressure;
                break;
            }
            pad->soundspeed = gdata->eos->getSoundSpeed(pad->pressure, density);
        }
    }

}

double PelletSolver::neon_radiation_power_density(double rho, double T){

    double temp,Tl,Tr,Z,Zl,Zr,L,Ll,Lr,k,n_rho;
    int j;
    // radiation power density = n^2 * Z * L * 1.e-15 [W/cm^3]

    n_rho = rho/3.350971e-23; // number density of neon nuclei

    //find lower index in the Neon radiation table
    temp = floor(10*T)/10.0;
    if (temp < 1.0)
        return 0.0;
    else if (temp < 5.0)
        j = round(10*(temp - 1));
    else if (temp <= 30.0)
        j = round((temp - 5)/0.5) + 40;
    else
        j = 90;

    Tl = neon_radiation_data(j,0);
    Tr = neon_radiation_data(j+1,0);

    Zl = neon_radiation_data(j,1);
    Zr = neon_radiation_data(j+1,1);

    Ll = neon_radiation_data(j,2);
    Lr = neon_radiation_data(j+1,2);

    k = (T - Tl)/(Tr - Tl);
    Z = (1-k)*Zl + k*Zr;
    L = (1-k)*Ll + k*Lr;

    return n_rho*n_rho*Z*L*1.e-17; // power density using ms time unit
}

double PelletSolver::neon_radiation_data(int i, int j){

    static  double rad[91][3];
    bool First = true;
    if(First){
        First = false;
        rad[0][0] =   1.000000e+00;
        rad[1][0] =   1.100000e+00;
        rad[2][0] =   1.200000e+00;
        rad[3][0] =   1.300000e+00;
        rad[4][0] =   1.400000e+00;
        rad[5][0] =   1.500000e+00;
        rad[6][0] =   1.600000e+00;
        rad[7][0] =   1.700000e+00;
        rad[8][0] =   1.800000e+00;
        rad[9][0] =   1.900000e+00;
        rad[10][0] =  2.000000e+00;
        rad[11][0] =  2.100000e+00;
        rad[12][0] =  2.200000e+00;
        rad[13][0] =  2.300000e+00;
        rad[14][0] =  2.400000e+00;
        rad[15][0] =  2.500000e+00;
        rad[16][0] =  2.600000e+00;
        rad[17][0] =  2.700000e+00;
        rad[18][0] =  2.800000e+00;
        rad[19][0] =  2.900000e+00;
        rad[20][0] =  3.000000e+00;
        rad[21][0] =  3.100000e+00;
        rad[22][0] =  3.200000e+00;
        rad[23][0] =  3.300000e+00;
        rad[24][0] =  3.400000e+00;
        rad[25][0] =  3.500000e+00;
        rad[26][0] =  3.600000e+00;
        rad[27][0] =  3.700000e+00;
        rad[28][0] =  3.800000e+00; 
        rad[29][0] =  3.900000e+00; 
        rad[30][0] =  4.000000e+00; 
        rad[31][0] =  4.100000e+00; 
        rad[32][0] =  4.200000e+00; 
        rad[33][0] =  4.300000e+00; 
        rad[34][0] =  4.400000e+00; 
        rad[35][0] =  4.500000e+00; 
        rad[36][0] =  4.600000e+00; 
        rad[37][0] =  4.700000e+00; 
        rad[38][0] =  4.800000e+00; 
        rad[39][0] =  4.900000e+00; 
        rad[40][0] =  5.000000e+00; 
        rad[41][0] =  5.500000e+00; 
        rad[42][0] =  6.000000e+00; 
        rad[43][0] =  6.500000e+00; 
        rad[44][0] =  7.000000e+00; 
        rad[45][0] =  7.500000e+00; 
        rad[46][0] =  8.000000e+00; 
        rad[47][0] =  8.500000e+00; 
        rad[48][0] =  9.000000e+00; 
        rad[49][0] =  9.500000e+00; 
        rad[50][0] =  1.000000e+01; 
        rad[51][0] =  1.050000e+01; 
        rad[52][0] =  1.100000e+01; 
        rad[53][0] =  1.150000e+01; 
        rad[54][0] =  1.200000e+01; 
        rad[55][0] =  1.250000e+01; 
        rad[56][0] =  1.300000e+01; 
        rad[57][0] =  1.350000e+01; 
        rad[58][0] =  1.400000e+01; 
        rad[59][0] =  1.450000e+01; 
        rad[60][0] =  1.500000e+01; 
        rad[61][0] =  1.550000e+01; 
        rad[62][0] =  1.600000e+01; 
        rad[63][0] =  1.650000e+01; 
        rad[64][0] =  1.700000e+01; 
        rad[65][0] =  1.750000e+01; 
        rad[66][0] =  1.800000e+01; 
        rad[67][0] =  1.850000e+01; 
        rad[68][0] =  1.900000e+01; 
        rad[69][0] =  1.950000e+01; 
        rad[70][0] =  2.000000e+01; 
        rad[71][0] =  2.050000e+01; 
        rad[72][0] =  2.100000e+01; 
        rad[73][0] =  2.150000e+01; 
        rad[74][0] =  2.200000e+01; 
        rad[75][0] =  2.250000e+01; 
        rad[76][0] =  2.300000e+01; 
        rad[77][0] =  2.350000e+01; 
        rad[78][0] =  2.400000e+01; 
        rad[79][0] =  2.450000e+01; 
        rad[80][0] =  2.500000e+01; 
        rad[81][0] =  2.550000e+01; 
        rad[82][0] =  2.600000e+01; 
        rad[83][0] =  2.650000e+01; 
        rad[84][0] =  2.700000e+01; 
        rad[85][0] =  2.750000e+01; 
        rad[86][0] =  2.800000e+01; 
        rad[87][0] =  2.850000e+01; 
        rad[88][0] =  2.900000e+01; 
        rad[89][0] =  2.950000e+01; 
        rad[90][0] =  3.000000e+01; 
        rad[0][1] =   2.855628e-02;
        rad[1][1] =   1.463662e-01;
        rad[2][1] =   4.330838e-01;
        rad[3][1] =   7.347827e-01;
        rad[4][1] =   8.955859e-01;
        rad[5][1] =   9.587141e-01;
        rad[6][1] =   9.825625e-01;
        rad[7][1] =   9.920914e-01;
        rad[8][1] =   9.962807e-01;
        rad[9][1] =   9.985018e-01;
        rad[10][1] =  1.000286e+00;
        rad[11][1] =  1.002707e+00;
        rad[12][1] =  1.006990e+00;
        rad[13][1] =  1.014886e+00;
        rad[14][1] =  1.028927e+00;
        rad[15][1] =  1.052489e+00;
        rad[16][1] =  1.089469e+00;
        rad[17][1] =  1.143327e+00;
        rad[18][1] =  1.215524e+00;
        rad[19][1] =  1.303975e+00;
        rad[20][1] =  1.402701e+00;
        rad[21][1] =  1.503367e+00;
        rad[22][1] =  1.597987e+00;
        rad[23][1] =  1.681105e+00;
        rad[24][1] =  1.750426e+00;
        rad[25][1] =  1.806171e+00;
        rad[26][1] =  1.849998e+00;
        rad[27][1] =  1.884091e+00;
        rad[28][1] =  1.910620e+00;
        rad[29][1] =  1.931495e+00;
        rad[30][1] =  1.948310e+00;
        rad[31][1] =  1.962368e+00;
        rad[32][1] =  1.974739e+00;
        rad[33][1] =  1.986323e+00;
        rad[34][1] =  1.997900e+00;
        rad[35][1] =  2.010166e+00;
        rad[36][1] =  2.023761e+00;
        rad[37][1] =  2.039274e+00;
        rad[38][1] =  2.057237e+00;
        rad[39][1] =  2.078115e+00;
        rad[40][1] =  2.102280e+00;
        rad[41][1] =  2.275239e+00;
        rad[42][1] =  2.498072e+00;
        rad[43][1] =  2.690690e+00;
        rad[44][1] =  2.819845e+00;
        rad[45][1] =  2.900079e+00;
        rad[46][1] =  2.954544e+00;
        rad[47][1] =  3.000644e+00;
        rad[48][1] =  3.049845e+00;
        rad[49][1] =  3.109150e+00;
        rad[50][1] =  3.181565e+00;
        rad[51][1] =  3.266204e+00;
        rad[52][1] =  3.358960e+00;
        rad[53][1] =  3.454218e+00;
        rad[54][1] =  3.546805e+00;
        rad[55][1] =  3.633317e+00;
        rad[56][1] =  3.712484e+00;
        rad[57][1] =  3.784763e+00;
        rad[58][1] =  3.851660e+00;
        rad[59][1] =  3.915086e+00;
        rad[60][1] =  3.976912e+00;
        rad[61][1] =  4.038694e+00;
        rad[62][1] =  4.101567e+00;
        rad[63][1] =  4.166216e+00;
        rad[64][1] =  4.232939e+00;
        rad[65][1] =  4.301736e+00;
        rad[66][1] =  4.372415e+00;
        rad[67][1] =  4.444695e+00;
        rad[68][1] =  4.518282e+00;
        rad[69][1] =  4.592917e+00;
        rad[70][1] =  4.668399e+00;
        rad[71][1] =  4.744582e+00;
        rad[72][1] =  4.821367e+00;
        rad[73][1] =  4.898686e+00;
        rad[74][1] =  4.976481e+00;
        rad[75][1] =  5.054694e+00;
        rad[76][1] =  5.133261e+00;
        rad[77][1] =  5.212099e+00;
        rad[78][1] =  5.291106e+00;
        rad[79][1] =  5.370152e+00;
        rad[80][1] =  5.449081e+00;
        rad[81][1] =  5.527706e+00;
        rad[82][1] =  5.605804e+00;
        rad[83][1] =  5.683127e+00;
        rad[84][1] =  5.759401e+00;
        rad[85][1] =  5.834343e+00;
        rad[86][1] =  5.907667e+00;
        rad[87][1] =  5.979103e+00;
        rad[88][1] =  6.048405e+00;
        rad[89][1] =  6.115366e+00;
        rad[90][1] =  6.179824e+00;
        rad[0][2] =   1.294421e-17;
        rad[1][2] =   5.032943e-17;
        rad[2][2] =   1.175009e-16;
        rad[3][2] =   1.618292e-16;
        rad[4][2] =   1.633444e-16;
        rad[5][2] =   1.487624e-16;
        rad[6][2] =   1.340048e-16;
        rad[7][2] =   1.251299e-16;
        rad[8][2] =   1.267496e-16;
        rad[9][2] =   1.459262e-16;
        rad[10][2] =  1.946532e-16;
        rad[11][2] =  2.921758e-16;
        rad[12][2] =  4.680815e-16;
        rad[13][2] =  7.667701e-16;
        rad[14][2] =  1.254301e-15;
        rad[15][2] =  2.028582e-15;
        rad[16][2] =  3.232404e-15;
        rad[17][2] =  5.064973e-15;
        rad[18][2] =  7.782106e-15;
        rad[19][2] =  1.167420e-14;
        rad[20][2] =  1.702276e-14;
        rad[21][2] =  2.405655e-14;
        rad[22][2] =  3.293476e-14;
        rad[23][2] =  4.376405e-14;
        rad[24][2] =  5.663246e-14;
        rad[25][2] =  7.163800e-14;
        rad[26][2] =  8.890362e-14;
        rad[27][2] =  1.085818e-13;
        rad[28][2] =  1.308546e-13;
        rad[29][2] =  1.559355e-13;
        rad[30][2] =  1.840726e-13;
        rad[31][2] =  2.155569e-13;
        rad[32][2] =  2.507330e-13;
        rad[33][2] =  2.900121e-13;
        rad[34][2] =  3.338860e-13;
        rad[35][2] =  3.829394e-13;
        rad[36][2] =  4.378594e-13;
        rad[37][2] =  4.993941e-13;
        rad[38][2] =  5.685169e-13;
        rad[39][2] =  6.461584e-13;
        rad[40][2] =  7.333081e-13;
        rad[41][2] =  1.340709e-12;
        rad[42][2] =  2.244150e-12;
        rad[43][2] =  3.337638e-12;
        rad[44][2] =  4.502229e-12;
        rad[45][2] =  5.694331e-12;
        rad[46][2] =  6.914699e-12;
        rad[47][2] =  8.175398e-12;
        rad[48][2] =  9.488725e-12;
        rad[49][2] =  1.086383e-11;
        rad[50][2] =  1.230390e-11;
        rad[51][2] =  1.380364e-11;
        rad[52][2] =  1.534886e-11;
        rad[53][2] =  1.691982e-11;
        rad[54][2] =  1.849633e-11;
        rad[55][2] =  2.006210e-11;
        rad[56][2] =  2.160670e-11;
        rad[57][2] =  2.312510e-11;
        rad[58][2] =  2.461607e-11;
        rad[59][2] =  2.608032e-11;
        rad[60][2] =  2.751923e-11;
        rad[61][2] =  2.893380e-11;
        rad[62][2] =  3.032439e-11;
        rad[63][2] =  3.169057e-11;
        rad[64][2] =  3.303126e-11;
        rad[65][2] =  3.434488e-11;
        rad[66][2] =  3.562955e-11;
        rad[67][2] =  3.688312e-11;
        rad[68][2] =  3.810312e-11;
        rad[69][2] =  3.928663e-11;
        rad[70][2] =  4.043005e-11;
        rad[71][2] =  4.152896e-11;
        rad[72][2] =  4.257794e-11;
        rad[73][2] =  4.357054e-11;
        rad[74][2] =  4.449942e-11;
        rad[75][2] =  4.535653e-11;
        rad[76][2] =  4.613346e-11;
        rad[77][2] =  4.682193e-11;
        rad[78][2] =  4.741422e-11;
        rad[79][2] =  4.790376e-11;
        rad[80][2] =  4.828558e-11;
        rad[81][2] =  4.855674e-11;
        rad[82][2] =  4.871657e-11;
        rad[83][2] =  4.876686e-11;
        rad[84][2] =  4.871173e-11;
        rad[85][2] =  4.855750e-11;
        rad[86][2] =  4.831223e-11;
        rad[87][2] =  4.798535e-11;
        rad[88][2] =  4.758706e-11;
        rad[89][2] =  4.712785e-11;
        rad[90][2] =  4.661804e-11;
    }
    return rad[i][j];
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