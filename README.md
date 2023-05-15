# 2D torus
Welcome to the readme. The goal of this code is to provide demos related to our paper [INSERT CITATION]
# Data
All of the files are stored as $n \times n \times 2$ matlab arrays of the velocity fields for snapshots, or as $n\times n\times2\times N_t$ for timeseries of length $N_t$. I.e. for a velocity field u, u(:,:,1) gives the x component of the field, and u(:,:,2) gives the y component. The length of time integration for converged Euler solutions can given by 4 times domain.Lt
# Demo codes
Currently, three demos are provided. 
* src/integrationDemo.m -- loads in a snapshot of turbulence saved on a $256\times256$ grid and integrates it on a $512\times512$ grid for 10 time units. Afterwards, it visualizes the trajectory.
* src/periodicOrbitConvergeDemo -- loads in a snapshot of turbulence saved on a $256\times256$, smoothes the initial condition, and then converges it as a pre-periodic orbit in Euler
* src/equilibriumConvergeDemo -- creates a good initial guess for converging an Euler equilibrium and converges it.
