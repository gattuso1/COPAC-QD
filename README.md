# QD-dimer

## Files for QD_quest

### Compile                   
Compiling command (./Compile), is gfortran but can be changed)
### Constants_au.f90             
Most constants used in the code (in atomic units) + conversion factors from a.u to s.i.
### Integrals.f90          
Contains all integral routines, RK, pulse routine ... 
### specfun.f90  
Contains some special functions
### Variables_au.f90
Define variables, allocatables, reads input .def and builds all system and pulses parameters
### Vectors.f90
Routines for randomly oriented vectors
### QD_quest.def
Input file
### QD_quest.f90
Main code, does everything, writes outputs, lauches dynamics, does matrix diagonalisation 
### Parameters-ex-55.f90
Contains parameters for 0.55 nm linker, exchange term
### Parameters-ex-02.f90
Contains parameters for 0.20 nm linker, exchange term
### Parameters-dir-55.f90
Contains parameters for 0.55 nm linker, direct term
### Parameters-dir-02.f90
Contains parameters for 0.20 nm linker, direct term
###Output.f90
File generating full outputs ... is not so good ... 
###Normal.f90
Generates good random numbers for normal distributions
###Make_Ham.f90
Creates all matrices (Hamiltonians, TDM)

## Files for computation of Coulomb integrals
### fetchIntegrals.f90      
### GetVarex.f90
### GetVar.f90  
### GetVecex.f90
### GetVec.f90  
### Integrals_Dexmax.f90   
### Integrals_Dmax.f90     

## Random files for other purposes
###Absorption.f90
###UnitPulse.f90
Converts pulse strength V/m to nJ

