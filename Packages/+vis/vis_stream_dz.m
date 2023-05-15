function vis_stream_dz(domain,state,tics,m)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% state  -- can either be a single velocity field or a trajectory itself
% pvis   -- percent vis. Lets say you want to view the trajectory at 20%
%           intervals, seeing 5 slices of the trajectory in total, 
%           then pvis = .2;


figure 

if size(state,4) > 1
    for j = 1:m:size(state,4)
        pause(0.0001);
    s = util.fftstream(domain,state(:,:,:,j));
ux = state(:,:,1,j);
uy = state(:,:,2,j);
colormap jet
% h = contourf(domain.X,domain.Y,s, 15,'EdgeColor','None');
imagesc(s);
hold on

xx = linspace(3,domain.Nx-2,tics);
yy = linspace(3,domain.Ny-2,tics);
[xx,yy] = meshgrid(xx,yy);
[X,Y] = meshgrid(1:domain.Nx,1:domain.Ny);
quiver(xx,yy,interp2(X,Y,ux,xx,yy),interp2(X,Y,uy,xx,yy),'k')
% h = pcolor(x,y,omega);
% set(h,'EdgeColor','none')
% hold on
% h = quiver(x(1:10:end,1:10:end),y(1:10:end,1:10:end),ux(1:10:end,1:10:end),uy(1:10:end,1:10:end));
% set(h,'Color','k')
colorbar
set(gca,'YDir','normal')
daspect([1,1,1]);

warning('on', 'MATLAB:contour:ConstantData')
    end
else

s = util.fftstream(domain,state);
ux = state(:,:,1);
uy = state(:,:,2);
colormap jet
% h = contourf(domain.X,domain.Y,s, 15,'EdgeColor','None');
imagesc(s);
hold on

xx = linspace(3,domain.Nx-2,tics);
yy = linspace(3,domain.Ny-2,tics);
[xx,yy] = meshgrid(xx,yy);
[X,Y] = meshgrid(1:domain.Nx,1:domain.Ny);
quiver(xx,yy,interp2(X,Y,ux,xx,yy),interp2(X,Y,uy,xx,yy),'k')
% h = pcolor(x,y,omega);
% set(h,'EdgeColor','none')
% hold on
% h = quiver(x(1:10:end,1:10:end),y(1:10:end,1:10:end),ux(1:10:end,1:10:end),uy(1:10:end,1:10:end));
% set(h,'Color','k')
colorbar
set(gca,'YDir','normal')
daspect([1,1,1]);

warning('on', 'MATLAB:contour:ConstantData')
end
end