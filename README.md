# Description
This repository contains sample code and data described in the paper "Exact Coherent Structures in Two-Dimensional Turbulence" by D. Zhigunov and R. O. Grigoriev.

# Data
The data files store snapshots of the flow field as $n \times n \times 2$ matlab arrays of the velocity fields. In particular, u(:,:,1) gives the x component of the velocity field, and u(:,:,2) gives the y component. Time series of snapshots are stored as $n\times n\times2\times N_t$ matlab arrays, where $N_t$ is the number of snapshots.

Other parameters are stored as elements of the domain object. Namely,
* domain.nu is the viscosity
* domain.hv is the array of parameters for the hyperviscos term: [k_v,alpha]
* domain.Lt is the time interval. Note that, for periodic orbits, this corresponds to a quarter of the temporal period, since these orbits are computed as preperiodic solutions with pi/2 rotation.
* domain.npu gives the number of timesteps per unit of time used in Newton-GMRES.

# Demo codes
Currently, three demos are provided.
* src/integrationDemo.m -- loads in a snapshot of turbulent flow saved on a $256\times256$ grid and integrates it on a $512\times512$ grid for 10 time units. Afterwards, the trajectory is visualized.
* src/periodicOrbitConvergeDemo -- loads in a snapshot of turbulent flow saved on a $256\times256$, applies smoothing to generate a corresponding initial condition for Newton-GMRES, and then converges it to a pre-periodic solution of the Euler equation.
* src/equilibriumConvergeDemo -- Creates an analytic initial guess that is close to an Euler equilibrium for Newton-GMRES and then converges it to an equilibrium of the Euler equation.
