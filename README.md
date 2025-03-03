Simulation of Cardiac Ion Currents

This repository contains Fortran source code for simulating cardiac ion currents. The model includes calculations for calcium and voltage-gated currents, along with stochastic binomial processes to model ion channel behavior.

Files in the Repository

binom.f - Subroutines for binomial random processes used in ion channel kinetics.

ca-currents.f - Implements calcium current dynamics.

v-currents.f - Implements voltage-gated current dynamics, including sodium currents.

eadca.f - The main program, which sets up and runs the simulation.


Prerequisites

To compile and run the code, you need the Intel Fortran Compiler (ifx). If it is not installed, you can get it from Intel's oneAPI toolkit:

Install Intel oneAPI

Compilation and Linking

To compile and link all files using ifx, navigate to the directory containing the source files and run the following command:

ifx eadca.f ca-currents.f v-currents.f binom.f -o eadca.out

Explanation of Compilation Steps:

ifx is the Intel Fortran Compiler.

eadca.f is the main simulation program.

ca-currents.f and v-currents.f contain calcium and voltage-gated ion channel dynamics.

binom.f provides helper functions for stochastic ion channel behavior.

-o eadca.out specifies the output executable name (eadca.out).

Running the Simulation

After compiling, run the executable:

./eadca.out

The program will execute the simulation and output the results based on the implemented model.

Customizing the Simulation

To modify the simulation parameters:

Grid Size: Defined in eadca.f with Lx and Ly.

Number of Stimuli: Controlled by the parameter nstim in eadca.f.

Ion Channel Properties: Adjust rates in ca-currents.f and v-currents.f.

License
