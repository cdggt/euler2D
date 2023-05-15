function vis_stream(domain,state,pvis)
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
        stream = util.fftstream(domain,state(:,:,:,j));
%         omega = util.fftsquaregridchange(omega,fluidparams.Nx,128);
       
        colormap parula
%         imagesc(stream);
        
        x = linspace(min(min(domain.X)),max(max(domain.X)),15);
y = linspace(min(min(domain.Y)),max(max(domain.Y)),15);
[x,y] = meshgrid(x,y);

        imagesc(x(1,:),y(:,1),stream)
        hold on
        quiver(x,y,interp2(domain.X,domain.Y,state(:,:,1,j),x,y),interp2(domain.X,domain.Y,state(:,:,2,j),x,y),'k')
        set(gca,'YDir','normal')
        axis equal
hold off;
        colorbar 
%         h = quiver(x(2:10:end-1,2:10:end-1),y(2:10:end-1,2:10:end-1),u(2:10:end-1,2:10:end-1),v(2:10:end-1,2:10:end-1),.5);
%         set(h,'Color','k')
        xlim(xbnd);
        ylim(ybnd);

        t=(j-1)*dt;
        title(['t = ' sprintf('%1.2f',t) '/'  sprintf('%1.2f',T), ' s']);

        pause(.25)
    end

else
    
figure 

stream = util.fftstream(domain,state);
ux = state(:,:,1);
uy = state(:,:,2);
colormap(BWR2())
h = contourf(domain.X,domain.Y,stream, 15,'EdgeColor','None');
% imagesc(stream);
hold on

xx = linspace(min(min(domain.X)),max(max(domain.X)),15);
yy = linspace(min(min(domain.Y)),max(max(domain.Y)),15);
[xx,yy] = meshgrid(xx,yy);
% [X,Y] = meshgrid(1:domain.Nx,1:domain.Ny);
quiver(xx,yy,interp2(domain.X,domain.Y,ux,xx,yy),interp2(domain.X,domain.Y,uy,xx,yy),'k')
% h = pcolor(x,y,omega);
% set(h,'EdgeColor','none')
% hold on
% h = quiver(x(1:10:end,1:10:end),y(1:10:end,1:10:end),ux(1:10:end,1:10:end),uy(1:10:end,1:10:end));
% set(h,'Color','k')
colorbar
xlim(xbnd);
ylim(ybnd);

warning('on', 'MATLAB:contour:ConstantData')
end

