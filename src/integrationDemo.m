clc; clear; close all; 
format long g
addpath('./../Packages')


load turbulenceSnapshot % load in initial condition
% set up domain
domain.Nx = 512; % gridsize in x
domain.Ny = 512; % gridsize in y
domain.nu = 10^-5; % viscosity
domain.Lt = 10; % length of time integration
domain.f = @(x,y) sin(4*x).*sin(4*y); % curl of the forcing

domain = dom.makeDomain(domain); % remake the domain for internal reasons
u = util.fftgridchange(u,size(u,1),domain.Nx); % change the size of u to match domain
u = util.dealias(domain,u,domain.Nx/3); % dealias u to match 2/3 dealiasing

t = linspace(0,domain.Lt,200); % define temporal grid for output
initStep = 1/2000;
minStep = 1/12345;
maxStep = 1/123;
tol = 1.7;
[u,hh] = int.rk45(domain,u,t,initStep,tol,minStep,maxStep);

% visualize the vorticity
vis.vis_omega(domain,u,1);