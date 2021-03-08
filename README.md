# Molecular Dynamics program
Simple molecular dynamics program to simulate particles interacting with a Lennard-Jones potential

Compile with the fortran2008 standard for support for the execute_command_line function call
This call allows the program to create an output directory

Code was written in 2016 during my postgraduate masters course and may not comply with the current Fortran standards

## Overview of the code 

1. Get the initial coordinates of the particles and other parameters
2. Initialise any other variables needed throughout the main program e.g loop counters, summation variables and variables for physical properties.
3. The main loop of the calculation which updates the positions and velocities of the particles. The positions and velocities of particles were updated after calculating the forces and then using those forces in the Velocity Verlet algorithm.
4. Output the data of interest e.g. the updated coordinates which could then be visualized in another program.


## Video of a stable cluster

[MD video](https://github.com/ianshepherd111/molecular_dynamics/blob/main/md.mp4)
