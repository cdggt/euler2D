function domain = makeDomain(options)

% Create parameters not present in the options
domain  = makeDefaults(options);
Nx = domain.Nx;
Ny = domain.Ny;

% generate grid data
x = (0:Nx-1)./Nx*2*pi; % define x vector corresponding to physical coordinates
y = (0:Ny-1)./Ny*2*pi; % define y vector corresponding to physical coordinates
[X, Y] = meshgrid(x,y);
domain.X = X;
domain.Y = Y;
domain.dx = 2*pi/Nx;
domain.dy = 2*pi/Ny;
domain.Nt = max(ceil(domain.Lt*domain.npu),4);

% define frequency grid
if(mod(Nx,2)==0)
    domain.qx = [0:(Nx/2)-1 (Nx/2) -(Nx/2)+1:-1]';
else
    domain.qx = [0:ceil(Nx/2)-1 -ceil(Nx/2)+1:-1]';
end
if(mod(Ny,2)==0)
    domain.qy = [0:(Ny/2)-1 (Ny/2) -(Ny/2)+1:-1]';
else
    domain.qy = [0:ceil(Ny/2)-1 -ceil(Ny/2)+1:-1]';
end

% Stored Parameters (Precalculated for speedup. These will be used often)
Dx = 1i*repmat(domain.qx',Ny,1);
Dy = 1i*repmat(domain.qy,1,Nx);
D2 = Dx.^2 + Dy.^2;	% Laplacian
Dx = fftshift(Dx);
Dy = fftshift(Dy);
Dy(1,:) = 0;
Dx(:,1) = 0;
Dx = ifftshift(Dx);
Dy = ifftshift(Dy);
ID2= 1./D2;  % invert the Laplacian (taking care of the central zero)
ID2(ID2==Inf) = -1e-32; % set the inverse's center to 0 rather than Inf
domain.stored = struct;
domain.stored.fftDx = Dx;
domain.stored.fftDy = Dy;
domain.stored.fftID2 = ID2;
domain.stored.fftD2 = D2;
domain.stored.fftcurlf = fft2(domain.amp*domain.f(domain.X,domain.Y));



end

function opts = makeDefaults(opts)

if ~isfield(opts,'N')
    opts.N = 512;
end

if ~isfield(opts,'Nx')
    opts.Nx = opts.N;
end

if ~isfield(opts,'Ny')
    opts.Ny = opts.N;
end

if ~isfield(opts,'Lt')
    opts.Lt = 1;
end

if ~isfield(opts,'npu')
    opts.npu = 300;
end

if ~isfield(opts,'nu')
    opts.nu = 1/100;
end


if ~isfield(opts,'f')
    opts.f = @(x,y) 0;
end

if ~isfield(opts,'A')
    opts.A = 1;
end

if ~isfield(opts,'hv')
    opts.hv = [0,0]; % wavenumber, slope
end

end