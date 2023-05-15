function wrap_omega(domain,state,pvis)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% state  -- can either be a single velocity field or a trajectory itself
% pvis   -- percent vis. Lets say you want to view the trajectory at 20%
%           intervals, seeing 5 slices of the trajectory in total, 
%           then pvis = .2;

percent = .25;
bufferx = round(domain.Nx*percent);
buffery = round(domain.Ny*percent);
% X = domain.X;
% Y = domain.Y;
% bufferLx = X(1,bufferx+1);
% bufferLy = Y(buffery+1,1);
% X = [X(:,1:bufferx)-bufferLx X X(:,end-bufferx+1:end)+bufferLx];
% X = [X(end-buffery+1:end,:); X; X(1:buffery,:)];
% Y = [Y(1:buffery,:)-bufferLy; Y; Y(end-buffery+1:end,:)+bufferLy];
% Y = [Y(:,end-bufferx+1:end) Y Y(:,1:bufferx)];

% xbnd = [X(1) X(end)];
% ybnd = [Y(1) Y(end)];

figure
N = size(state,4);
T = domain.Lt;
dt = T/(N-1);

if(N>1)
    mov = true;
    if(~exist('pvis','var'))
        pvis = 1/N;
    end
    visidx = round(linspace(1,N,round(1/pvis)));
    M = numel(visidx);
else
    mov = false;
    M = 1;
    visidx = [1];
end

    
for i = 1:M
    
    clf
    j = visidx(i);
    
    omega = util.fftcurl(domain,state(:,:,:,j));
    omega = [omega(:,end-bufferx+1:end) omega omega(:,1:bufferx)];
    omega = [omega(end-buffery+1:end,:); omega; omega(1:buffery,:)];
    
    contourf(omega, 15,'EdgeColor','None');
    hold on
    plot([1 domain.Nx+2*bufferx],[buffery buffery],'--k');
    plot([1 domain.Nx+2*bufferx],[domain.Ny+buffery domain.Ny+buffery],'--k');
    plot([bufferx bufferx],[0 domain.Ny+2*buffery],'--k');
    plot([domain.Nx+bufferx domain.Nx+bufferx],[0 domain.Ny+2*buffery],'--k');
    
    if(i==1 || mov)
        colormap parula
        colorbar
        xlim([0 domain.Nx+2*bufferx]);
        ylim([0 domain.Ny+2*buffery]);
        set(gca,'xticklabel',{[]})
        set(gca,'yticklabel',{[]})
    end
    
    if(mov)
        t=(j-1)*dt;
        title(['t = ' sprintf('%1.2f',t) '/'  sprintf('%1.2f',T), ' s']);
    end
    
    pause(.25)
    
end
    
warning('on', 'MATLAB:contour:ConstantData')

end