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
## Data Structure: 
STL Vector.  
Consider testing speed/memory using vector.  
The std::vector support random access which allows constant-time access to any element.  
Could use `reserve` or `resize` to manage memory allocation for vectors.

## Usage:

## Current Work:
We are currently working on the following jobs:  

**1. Particle data structure: James**

**2. Boundary Conditions: James**

**3. Laxwendroff solver: Chao**  
Working on the compute spatial derivatice and time integration.  
**4. 1D solver: Chao**  
Finished the code structure.  
**4. material library: Chao**  
Finished the basic neon class, get(1+Z*) need to be finished

## Future Work:
1. heat decomposition with heating source  
2. output result
3. use of input file
4. restart
5. makefile
6. time anaylsis
7. material 

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

