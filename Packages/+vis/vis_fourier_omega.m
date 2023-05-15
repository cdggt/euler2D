function vis_fourier_omega(domain,state,pvis)
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
       
        
        imagesc((abs(fftshift(fft2(omega)))));
%         colormap bone
%         imagesc(x(1,:),y(:,1),omega)
%         hold on
%         axis equal
        colorbar 
%         h = quiver(x(2:10:end-1,2:10:end-1),y(2:10:end-1,2:10:end-1),u(2:10:end-1,2:10:end-1),v(2:10:end-1,2:10:end-1),.5);
%         set(h,'Color','k')
        colormap parula
        t=(j-1)*dt;
        title(['t = ' sprintf('%1.2f',t) '/'  sprintf('%1.2f',T), ' s']);

        pause(.1)
    end

else
    
figure 

omega = util.fftcurl(domain,state);
% colormap bone
imagesc(abs(fftshift(fft2(omega))));
% h = pcolor(x,y,omega);
% set(h,'EdgeColor','none')
% hold on
% h = quiver(x(1:10:end,1:10:end),y(1:10:end,1:10:end),ux(1:10:end,1:10:end),uy(1:10:end,1:10:end));
% set(h,'Color','k')
colorbar

warning('on', 'MATLAB:contour:ConstantData')
end

