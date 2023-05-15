# 2D torus
Welcome to the readme. The goal of this code is to provide demos related to our paper [INSERT CITATION]
# Data
For single snapshots in time, all of the files are stored as $n \times n \times 2$ velocity fields. I.e. for a velocity field u, u(:,:,1) gives the x component of the field, and u(:,:,2) gives the y component. Trajectories of both turbulence and Euler, as output by DNS, are given as $n\times n\times2\times N_t$ arrays for timeseries of length $N_t$.

Other parameters, such as viscosity, hyperviscous parameters, and length of time integration are given inside of the domain object. Namely,
* domain.nu is the viscosity
* domain.hv gives the array for hyperviscosity: [k_v,alpha]
* domain.Lt gives the length of time integration. Note that for periodic orbits this corresponds to a quarter of the period, since they are computed as preperiodic with pi/4 rotation.
* domain.npu gives the number of timesteps per unit used in Newton-GMRES.
# Demo codes
Currently, three demos are provided. 
* src/integrationDemo.m -- loads in a snapshot of turbulence saved on a $256\times256$ grid and integrates it on a $512\times512$ grid for 10 time units. Afterwards, it visualizes the trajectory.
* src/periodicOrbitConvergeDemo -- loads in a snapshot of turbulence saved on a $256\times256$, smoothes the initial condition, and then converges it as a pre-periodic orbit in Euler.
* src/equilibriumConvergeDemo -- creates a good initial guess for converging an Euler equilibrium and converges it.
