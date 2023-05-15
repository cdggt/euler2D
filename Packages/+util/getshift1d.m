function [shift] = getshift1d(u,v)

% if v = u(x-L), this finds L

u = u/norm(u);
v = v/norm(v);

U = fft(u);
V = fft(v);

F = fft2(conj(V).*U);

percentdistance = max(F(:));
shift = find(F==percentdistance,1); %this returns the INDEX of the shift

% EXAMPLE
% x = linspace(0,2*pi,101);
% x = x(1:end-1);
% signal1 = sin(x);
% signal2 = sin(x-.2513);
% index = getshift1d(signal1,signal2);
% L = x(index); % equals .2513

end