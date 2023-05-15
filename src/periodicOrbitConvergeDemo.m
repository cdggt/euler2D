clc; clear; close all; 
format long g
addpath('./../Packages')

%
load turbulenceSnapshot

% set up domain
domain.Nx = 256; % gridsize in x
domain.Ny = 256; % gridsize in y
domain.nu = 0; % viscosity
domain.f = @(x,y) 0; % curl of the forcing
domain = dom.makeDomain(domain);

u = util.fftgridchange(u,size(u,1),domain.Nx); % change grid of initial condition
u = util.dealias(domain,u,256/4); % lightly dealias initial condition
u = util.streamFunctionSmooth(domain,0.5,u); % apply stream function smoothing
u = util.fourierMask(domain,u,@(k) exp(-(k/7).^2)); % gaussian blur

% hyperviscous blurring
domain.hv = [10,0.5]; % set hyperviscosity with [k_v,alpha]
u = u*sqrt(5/util.Energy(u)); % renormalize energy
t = linspace(0,26,2);
[p,~] = int.rk45(domain,u,t,1/2000,1.7,1/12345,1/123);
u = p(:,:,:,end);

% re-renormalize
u = u*sqrt(5/util.Energy(u));

% find period of periodic orbit
domain.hv = [64,0.5];
t = linspace(0,26,26*120);
[p,hh] = int.rk45(domain,u,t,1/2000,1.7,1/12345,1/123);
d = zeros(size(t));
for j = 1:size(p,4)
    d(j) = sum(sum(sum(abs(util.fftcurl(domain,p(:,:,:,1)-p(:,:,:,j))))));
end
clear('p')
d(t<5) = inf;
[~,j] = min(d);
domain.Lt = t(j)/4; % since these are preperiodic orbits, we use T/4 for newton
domain.npu = round(mean(1./hh))+33; % number of points per timestep

% shift initial condition, which is necessary for rotation in Newton
[xPhase,yPhase] = util.findPhase(domain,u); % finds phase [-phi_x,-phi_y]
u = sym.transx(domain,sym.transy(domain,u,yPhase),xPhase);
% converge
[u,~,n] = newtonmethods.convergePPO(domain,u,'temp.mat',5*10^-5,101);

% integrate and visualize
domain.Lt = domain.Lt*4;
t = linspace(0,domain.Lt,round(domain.Lt*11));
[p,hh] = int.rk45(domain,u,t,1/2000,1.7,1/12345,1/123);
vis.vis_omega(domain,p,1);

