clc; clear; close all
format long g
addpath('./../Packages')

% intiialize domain
domain.nu = 0*10^-5;
domain.amp = 0*1;
domain.f = @(x,y) sin(4*x).*sin(4.*y);
domain.Lt = 0.1;
domain.Nx = 256;
domain.Ny = 256;
domain.npu = 222;
domain.hv = [0,0];
domain = dom.makeDomain(domain);

% we're going to artificially construct an initial condition here
x = domain.X-pi; % extract x and y coordinates
y = domain.Y-pi;
% construct vorticity field
w = (sin(x)+sin(y));
w = sign(w).*abs(w).^3.5;
w = util.dealias(domain,w,64); % dealias to remove potential high frequencies
u = util.fftuncurl(domain,w);
u = u*sqrt(5/util.Energy(u));
u = newtonmethods.convergeEq(domain,u,'temp.mat',10^-5);