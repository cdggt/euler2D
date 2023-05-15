function vis_omega(domain,p,k)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% state  -- can either be a single velocity field or a trajectory itself


x = domain.X;
y = domain.Y;
xbnd = [x(1) x(end)];
ybnd = [y(1) y(end)];


if(size(p,4)>1)
    figure()
    t = linspace(0,domain.Lt,size(p,4));
    p = p(:,:,:,1:k:end);
    t = t(1:k:end);
    for j = 1:length(t)
        u = util.fftgridchange(p(:,:,:,j),size(p,1),domain.Nx);
    w = util.fftcurl(domain,u);

    imagesc(w,[-max(max(abs(w))),max(max(abs(w)))]); hold on
    colormap(BWR2())
    title(num2str(t(j)));
%     colorbar
    set(gca,'YDir','normal')
        set(gca,'XTick',[])
    set(gca,'YTick',[])
%      contour(X,Y,w,[-10,0,10],'k'); hold off
    daspect([1,1,1]); hold off
    pause(0.0001);
    end

else
    
figure 
p = util.fftgridchange(p,size(p,1),domain.Nx);
w = util.fftcurl(domain,p);

    imagesc(w,[-max(max(abs(w))),max(max(abs(w)))]); hold on
    colormap(BWR2())
%     colorbar
    set(gca,'YDir','normal')
        set(gca,'XTick',[])
    set(gca,'YTick',[])
%      contour(X,Y,w,[-10,0,10],'k'); hold off
    daspect([1,1,1]); hold off

warning('on', 'MATLAB:contour:ConstantData')
end

