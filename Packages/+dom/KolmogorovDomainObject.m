function domain = KolmogorovDomainObject(varargin)

% Create parser
params = inputParser;
params = dom.default_kolmogorov_scheme(params);
params.parse(varargin{:});

domain = params.Results;

Lx = domain.Lx;  % width of the domain we are integrating over
Ly = domain.Ly;  % height of the domain we are inegrating over
Nx = domain.Nx;  % number of grid points x
Ny = domain.Ny;  % number of grid points y
Kx = domain.Kx;  % x forcing wave number 
Ky = domain.Ky;  % y forcing wave number (only used if setup is 'C')

%-------------------------------------------------------------------------
x = (0:Nx-1)./Nx.*Lx; % define x vector corresponding to physical coordinates
y = (0:Ny-1)./Ny.*Ly; % define y vector corresponding to physical coordinates
[X, Y] = meshgrid(x,y);

% forcing profiles in functional form (no amplitudes)
if domain.setup == 'K'
    f(:,:,1) = sin((2*pi*Ky/Ly).*Y);
    f(:,:,2) = 0.*f(:,:,1);
    curlf = -(2*pi*Ky/Ly).*cos((2*pi*Ky/Ly).*Y); % if f=sin(n y)xhat, then curl.f =-n cos(ny)zhat
else
    warning('Setup not supported')
    return;
end

% Domain Parameters
% domain.f = f;
domain.X = X;
domain.Y = Y;
domain.dx = Lx/Nx;
domain.dy = Ly/Ny;
domain.Nt = max(ceil(domain.Lt*domain.npu),4);

if(mod(Nx,2)==0)
    domain.qx = [0:(Nx/2)-1 (Nx/2) -(Nx/2)+1:-1]'.*2*pi/Lx;
else
    domain.qx = [0:ceil(Nx/2)-1 -ceil(Nx/2)+1:-1]'.*2*pi/Lx;
end
if(mod(Ny,2)==0)
    domain.qy = [0:(Ny/2)-1 (Ny/2) -(Ny/2)+1:-1]'.*2*pi/Ly;
else
    domain.qy = [0:ceil(Ny/2)-1 -ceil(Ny/2)+1:-1]'.*2*pi/Ly;
end

% Stored Parameters (Precalculated for speedup. These will be used often)
if(domain.kmax <= 0), domain.kmax = inf; end
Dx = 1i*repmat(domain.qx',Ny,1);
Dy = 1i*repmat(domain.qy,1,Nx);
D2 = Dx.^2 + Dy.^2;	% Laplacian
Dx = fftshift(Dx);
Dy = fftshift(Dy);
Dy(1,:) = 0;
Dx(:,1) = 0;
Dx = ifftshift(Dx);
Dy = ifftshift(Dy);
Dx(sqrt(abs(D2))>(domain.kmax))=0;
Dy(sqrt(abs(D2))>(domain.kmax))=0;
D2(sqrt(abs(D2))>(domain.kmax))=0;
ID2= 1./D2;  % invert the Laplacian (taking care of the central zero)
ID2(ID2==Inf) = -1e-32; % set the inverse's center to 0 rather than Inf
domain.stored.fftDx = Dx;
domain.stored.fftDy = Dy;
domain.stored.fftID2 = ID2;
domain.stored.fftD2 = D2;
domain.stored.fftcurlf = fft2(curlf);

% Physical parameters to compare to experiment
domain.physical = struct;
domain.physical.nmx = max(Kx*2,1); % Number of magnets in the x direction
domain.physical.nmy = max(Ky*2,1); % Number of magnets in the y direction
domain.physical.ngpx = Nx/domain.physical.nmx; % Number of grid points per magnet in the x direction
domain.physical.ngpy = Ny/domain.physical.nmy; % Number of grid points per magnet in the y direction
domain.physical.aspectratio = Ly/Lx;
domain.physical.Re = sqrt(domain.amp)/domain.nu*(Ly/(2*pi))^(2/3);

domain.append = @(varargin) append(domain,varargin{:}); %%% check to make sure this doesnt append params always to the original 
end

function domain = append(domain,varargin)
    nf = numel(varargin);
    if(nf == 1)
        params = varargin{1};
    else
        params = struct;
        for i = 1:2:numel(varargin)
            params.(varargin{i}) = varargin{i+1};
        end
    end
    fields = fieldnames(params);
    for i=1:numel(fields)
        field = fields{i};
        value = params.(field);
        switch(lower(field))
            case 'lt'
                domain.Lt = max(0,value);
                domain.Nt = ceil(domain.Lt*domain.npu);
            case 'phase'
                domain.phase = value;
            case 'nu'
                domain.nu = value;
            case 'beta'
                domain.beta = value;
            case 'alpha'
                domain.alpha = value;
            case 'npu'
                domain.npu = max(1,value);
                domain.Nt = ceil(domain.Lt*domain.npu);
            otherwise
                warning('field %s was not recognized and not appended',field);
        end
    end
    domain.append = @(varargin) append(domain,varargin{:});
end
