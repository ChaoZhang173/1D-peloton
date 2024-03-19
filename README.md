# 1D-peloton

This repository is used for developing 1D peloton code.  
The 1D peloton code is used do quick analysis for pellet ablation in tokamaks.

## Author:

Chao Zhang   
James Corbett

## Algorithm:
1. MHD equations:  
   Now using lax-wendroff and upwind (only for right boundar) scheme.    
   Other schemes maybe added in the future.
   The upwind is applied on the right boundary, laxwendroff has high oscillation near the boundary   
3. spatial derivative:  
   Use Newton Interpolation  
   unstable, may fail if some particles too close  
4. eos:
   polytropic eos, will add saha equations in the future
5. FIST:
   Currently use a ifStart to indicate if at the first step  
6. write result:
   (to be done)write initial state at the beginning while write the output function   
   other steps will be written normally.
7. local spacing:
   the local spaing = 0.5((distance to left) + (distance to right))  
   1st particle: local spacing = x_2nd - x_1st  
   last particle: local spacing = x_last - x_(last-1)  
9. electron density integral and heat deposition:   
   The 1D code only has left integral, and right integral is not used.  
   at the end of computeDensityIntegral, set the `-∇·q` and `q+-` to 0 for every partilce,  
   This is for the following function: computeHeatDeposition.
   1+Z* is related to Tinf, need to take care of in setMaterial fucntion  
   Currently the two functions are seperately, but they could merge together.  
   Maybe work on it later  
10. the state library is used to initializing the particles at the beginning,  
   in the following steps new particles will be added by generateBoundaryParticles in pellet_solver  
11. Generate boundary particles:  
   ~Currently generate particles at 0 - smallest x,~  
    Now use 0 - pelletV*dt  
    mass_fix = dx *Vinflow(100), particles have same mass  
    Note that the smallest x maybe incorrect, because didn't reorder the particles at last step  
13. dq/dx (QplusminusGradient):  
    dq/dx = (q_right - q_left)/(x_right - x_left)   
14. mindx:  
    The mindx in Global_Data is the smallest distance between two particles  
    If localspacing < 3*mindx will report a warning message and put particle to the middle
15. Right Boundary:
    The lax-wendroff method will introduce oscillation on the right boundary
    So change to use upwind method.  
       
## Data Structure: 
1. STL Vector.  
-  Consider testing speed/memory using vector.  
   The std::vector support random access which allows constant-time access to any element.  
   Could use `reserve` or `resize` to manage memory allocation for vectors.  
-  Switch to use smart pointer `std::unique_ptr<std::vector<pdata_t>> particle_data`  
   Smart pointer could manage memory automatically.  
2. Particle Data
-  Use a smart pointer particle_data to store all particles in Global_Data class.
-  The 1st particle is near but NOT at the pellet surface  
   Pellet surface: -0.5*initial spacing  
   New particles will be push and added into the particle list  
-  Neighbour: Left neighbour for p0 and right neighbour for plast   
   Local spaicng: 0.5((distance to left) + (distance to right))
3. Ghost Particle
-  A smart pointer to store all (currently 2) ghost particles in Global_Data
-  1st element(0):
   The left neighbour of 1st particle: x = x_1st - localspacing_1st;  
   parameter = 1st particle  
   It is symmetric about pellet (x_1st - 0.5*localspacing_1st)
-  2nd element(1):
   The right neighbour of last particle: x = x_last + 0.5*localspacing_last;
   parameter = vacumm  

## Usage:
-  The 1D peloton code uses C++20 for smart pointer (added at C++17)
-  Before compile, change initial seetings in the `initializer.cpp`  
-  compile by
```
make
```
delete by 
```
make clean
```  
-  run the code by
```
./lp <output_filename>
``` 
## Current Work:
We are currently working on the following jobs:  

**0. debugging: Chao**  
  Finished debugging, code could run successfully 03/19/2024  
  Fixed issue: uninitialized variable: vector<double> neeff 03/19/2024  
  The accuracy needs to be verified   
**1. timer: Chao**   

## Future Work:
1. use of input file
2. restart
3. merge the computeDensityIntegral with computeHeatDeposition 
4. release memory  

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
- geometry_pellet(not used)  
- boundary_pellet
- eos
- material_lib

2. EOS class: Chao 02/29/2024  
- Finished the EOS file:  
- The polytropic gas eos with EOSChoice = 1.

3. Laxwendroff solver: Chao 03/04/2024
- Finished spatial derivative.
- Finished the time integration.
4. 1D solver: Chao  03/04/2024    
- Finished the code structure.
- Finished the updatePatcileState.
- Finished moveParticle().
- Finished computeTemperature().
- Finished adjustDtByWriteTimeInterval().  
5. heating model: Chao 03/05/2024
- finished the computeDensityIntegral
- finished addHeatSource for heat deposition
- finished computeHeatDeposition  
6. Generate Boundary particles: Chao 03/11/2024    
- finished the initializing particles  
- finished genreateBoundaryParticle:03/11/2024  
- finished coumpteMassFlowRate  
Need to understand Pinflow and Vinflow
7. Particle data structure: Chao 03/11/2024  
- particle data/neighbour particle use smart pointer  
8. Material library: Chao  03/12/2024  
9. Initializer: Chao 02/12/2024     
- Initializer is partially finished with using user set input.  
- Finished the basic neon class   
- The 1+Z* function return 4.9
10. Boundary Condition: Chao 03/12/2024    
- finished compute boundary condition
- finished computeQplusminusGradient
11. Ghost and neighbour particles: Chao 03/14/2024  
- Finished generating ghost particles  
- Finished set neighbours
12. Reorder and re locate particles, update local spacing: Chao 03/14/2024  
- Finished reorder particles  
- Finished update local spacing  
- Finished re locate particles in updatelocalSpacing
13. Radiation cooling: Chao 03/14/2024
- Finished neon radiation cooling, used for eoschoice = 2 (saha neon)  
  but currently not used for eoschoice = 1 (polytropic)  
14. Set up eos and pellet material: Chao 03/15/2024  
- Finished pellet material   
- Finished setEOS for polytropic eos,   
  other eos will be added in the future
15. Output result: Chao 03/15/2024  
- Finished output result
- Write result after adjustDtByWriteTimeInterval() function
  Reason: don't need to reorder particles and update localspacing twice after moveParticle  
  However, the smallest x output by generateBoundaryParticle maybe incorrect  
  It will write a result at the 1st time step  
  By setting currenttime <= tend (instead of "<"), it will output at the result at last+1 step
16. Makefile: Chao 03/15/2024
- Finished makefile     
17. upwind solver for right boundary: Chao 03/19/2024  
- Finished upwind solver  
