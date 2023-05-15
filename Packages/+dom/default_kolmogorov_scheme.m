function [params] = default_kolmogorov_scheme(params)

%forcing parameters
params.addParamValue('amp', 1.0, @isnumeric); % Amplitude of forcing
params.addParamValue('Kx', 0.0, @isnumeric); % forcing wavenumber
params.addParamValue('Ky', 4.0, @isnumeric); % forcing wavenumber

%spatial domain parameters
default_Nx = 128;
default_Ny = 128;
default_Lx = 2*pi;
default_Ly = 2*pi;
% default_qx = [0:default_Nx/2-1 0 -default_Nx/2+1:-1]'.*2*pi/default_Lx;
% default_qy = [0:default_Ny/2-1 0 -default_Ny/2+1:-1]'.*2*pi/default_Ly;
params.addParamValue('setup', 'K', @(x) ismember(x,{'K','C'})); % type of forcing pattern
params.addParamValue('Lx', default_Lx, @(x) (x>0)); % physical length of domain in x direction
params.addParamValue('Ly', default_Ly, @(x) (x>0)); % physical length of domain in y direction
params.addParamValue('Nx', default_Nx, @(x) (x>0)&&(x-round(x)==0)); % number of grid points in x direction
params.addParamValue('Ny', default_Ny, @(x) (x>0)&&(x-round(x)==0)); % number of grid points in y direction
params.addParamValue('kmax', Inf, @isnumeric); % number of grid points in y direction

%temporal domain parameters
params.addParamValue('npu', 2^6, @(x) (x>0)&&(x-round(x)==0)); %number of time steps per unit of non-dimensional time
params.addParamValue('Lt',1, @(x) x>=0); % default integration time

%fluid parameters
params.addParamValue('nu', 0.0250, @isnumeric); %effective kinematic viscosity
params.addParamValue('alpha', 0.0, @isnumeric); %Rayleigh friction coefficient
params.addParamValue('beta', 1.0, @isnumeric);  %Advection Coefficient

%symmetry parameters
params.addParamValue('phase', 0.0, @isnumeric);  %phase shift along conintuous translational symmetry

end

