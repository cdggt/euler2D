function [shiftx,shifty] = getshift2d(u,v)

% if v = u(x-L), this finds L

u = u/norm(u(:));
v = v/norm(v(:));

U = fft2(u);
V = fft2(v);

F = fft2(conj(V).*U);

percentdistance = max(F(:));
[shifty,shiftx] = find(F==percentdistance,1); %this returns the INDEX of the shift

% EXAMPLE
% x = linspace(0,2*pi,101);
% x = x(1:end-1);
% [X,Y] = meshgrid(x,x);
% signal1 = sin(X).*cos(Y);
% signal2 = sin(X-.25).*cos(Y-.125);
% [shiftx,shifty] = getshift2d(signal1,signal2);
% Lx = x(shiftx) % equals .25
% Ly = x(shifty) % equals .125

end