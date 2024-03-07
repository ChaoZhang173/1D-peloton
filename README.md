# 1D-peloton

This repository is used for developing 1D peloton code.  
The 1D peloton code is used do quick analysis for pellet ablation in tokamaks.

## Author:

Chao Zhang   
James Corbett

## Algorithm:
1. MHD equations:  
   Now only using lax-wendroff scheme.  
   Other schemes maybe added in the future.
2. spatial derivative:  
   Use Newton Interpolation
3. eos:
   polytropic eos, will add saha equations in the future
4. FIST:
5. write result:
   (to be done)write initial state at the beginning while write the output function   
   other steps will be written normally.
6. local spacing:
   the local spaing = 0.5((distance to left) + (distance to right))
   no local spacing for the boundary particles
7. electron density integral and heat deposition:   
   the 1D code only has left integral, and right integral is not used.  
   at the end of computeDensityIntegral, set the `-∇·q` and `q+-` to 0 for every partilce,  
   this is for the following function: computeHeatDeposition.  
   Currently the two functions are seperately, but they could merge together.
   Maybe work on it later
8. the state library is used to initializing the particles at the beginning,  
   in the following steps new particles will be added by generateBoundaryParticles in pellet_solver  
   
       
## Data Structure: 
STL Vector.  
Consider testing speed/memory using vector.  
The std::vector support random access which allows constant-time access to any element.  
Could use `reserve` or `resize` to manage memory allocation for vectors.  
Switch to use smart pointer `std::unique_ptr<std::vector<pdata_t>> particle_data`  
smart pointer could manage memory automatically.  

## Usage:

## Current Work:
We are currently working on the following jobs:  

**1. Particle data structure: Chao**  
use smart pointer  
**2. Boundary Conditions: Chao**    
finished the initializing particles
**3. material library: Chao**  
     Finished the basic neon class, get(1+Z*) need to be finished   
     The 1+Z* function in 3D code returns 0, need to work with   
**4. Initializer: Chao**   
   Initializer is partially finished with using user set input.
   
## Future Work:
0. initialize the variables (important)
1. makefile
2. output result
3. use of input file
4. restart
5. time anaylsis
6. merge the computeDensityIntegral with computeHeatDeposition
7. set up eos
8. understand how state works

## Previous Work:
1. Code Structure: Chao  02/28/2024  
Finished the basic code structure (with content empty):
- lp_main
- initializer  
- particle_data   
- lp_solver   
- particle_viewer  
- pellet_solver  
- state_pellet  
- geometry_pellet  
- boundary_pellet
- eos
- material_lib

2. EOS class: Chao 02/29/2024  
Finished the EOS file:  
The polytropic gas eos with EOSChoice = 1.

3. Laxwendroff solver: Chao 03/04/2024  
     Finished spatial derivative.
     Finished the time integration.
4. 1D solver: Chao  03/04/2024    
     Finished the code structure.  
     Finished the updatePatcileState.  
     Finished moveParticle().  
     Finished computeTemperature().  
     Finished adjustDtByWriteTimeInterval().  
5. heating model: Chao 03/05/2024
- finished the computeDensityIntegral
- finished addHeatSource for heat deposition
- finished computeHeatDeposition
