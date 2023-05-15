function [xPhase,yPhase] = findPhase(domain,p)
options = struct;
options.Lx = 2*pi;  % width of the domain we are integrating over
options.Ly = 2*pi;  % height of the domain we are inegrating over
options.Ky = 2;     % this does nothing
options.Nx = size(p,2);   % spatial grid points along x
options.Ny = size(p,1);   % spatial grid points along y
options.nu = 10.^(-5); % viscosity 10^-4 -- 10^-6
options.amp = 1; % forcing amplitude -- make it whatever
options.Lt = 188;
options.npu = 1;
domain = dom.KolmogorovDomainObject(options);

xPhase = zeros(1,size(p,4));
yPhase = zeros(1,size(p,4));
for j = 1:size(p,4)
    s = util.fftstream(domain,p(:,:,:,j));
%     A1 = sum(sum(s.*cos(x)))*4*pi*pi/numel(s);
%     A2 = sum(sum(s.*cos(y)))*4*pi*pi/numel(s);
%     B1 = sum(sum(s.*sin(x)))*4*pi*pi/numel(s);
%     B2 = sum(sum(s.*sin(y)))*4*pi*pi/numel(s);
    S = fft2(s);
    xPhase(j) = angle(S(1,2));
    yPhase(j) = angle(S(2,1));
    
end
    
end