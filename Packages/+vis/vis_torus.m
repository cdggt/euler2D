function vis_torus(domain,state,pvis)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% state  -- can either be a single velocity field or a trajectory itself
% pvis   -- percent vis. Lets say you want to view the trajectory at 20%
%           intervals, seeing 5 slices of the trajectory in total, 
%           then pvis = .2;

x = domain.X;
y = domain.Y;
xbnd = [x(1) x(end)];
ybnd = [y(1) y(end)];

if(size(state,4)>1)
    
    om_levels = 20;
    N = size(state,4);
    if(~exist('pvis','var'))
        pvis = 1/N;
    end
    visidx = round(linspace(1,N,round(1/pvis)));
    M = numel(visidx);

    figure
    N = size(state,4);
    T = domain.Lt;
    dt = T/(N-1);
    for i = 1:M
        clf
        j = visidx(i);
        
%         u = state(:,:,1,i);
%         v = state(:,:,2,i);
        omega = util.fftcurl(domain,state(:,:,:,j));
%         omega = util.fftsquaregridchange(omega,fluidparams.Nx,128);
       
        R=1;
        r=.5;
        u=linspace(0,2*pi,domain.Nx);
        v=linspace(0,2*pi,domain.Ny);
        [u,v]=meshgrid(u,v);
        x=(R+r.*cos(v)).*cos(u); 
        y=(R+r.*cos(v)).*sin(u);
        z=r.*sin(v);
%         figure(1)
        surf(x,y,z,omega,'Edgecolor','none');
        % view([-52,64])
        colormap parula
        axis equal
        colorbar
        t=(j-1)*dt;
        title(['t = ' sprintf('%1.2f',t) '/'  sprintf('%1.2f',T), ' s']);

        pause(1)
    end

else
    
figure 

R=1;
r=.5;
u=linspace(0,2*pi,domain.Nx);
v=linspace(0,2*pi,domain.Ny);
[u,v]=meshgrid(u,v);
x=(R+r.*cos(v)).*cos(u); 
y=(R+r.*cos(v)).*sin(u);
z=r.*sin(v);
figure(1)
surf(x,y,z,util.fftcurl(domain,state),'Edgecolor','none');
% view([-52,64])
colormap parula
axis equal
colorbar

warning('on', 'MATLAB:contour:ConstantData')
end