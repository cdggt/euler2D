function vis_omega2000(domain,u)
warning('off', 'MATLAB:contour:ConstantData')
% domain -- domain object that the data was created on
% state  -- can either be a single velocity field or a trajectory itself


domain.Nx = max(size(u,1),2000);
domain.Ny = max(size(u,1),2000);
domain = dom.makeDomain(domain);


    
% figure('Units', 'pixels','outerposition',[600 400 340 400])
u = util.fftgridchange(u,size(u,1),domain.Nx);
w = util.fftcurl(domain,u);
x = linspace(0,2*pi,domain.Nx+1);
x = x(1:end-1);
x = x-pi;
    imagesc(x,x,w,[-max(max(abs(w))),max(max(abs(w)))]); hold on
    colormap(BWR2())
%     colorbar
    set(gca,'YDir','normal')
        set(gca,'XTick',[])
    set(gca,'YTick',[])
%      contour(X,Y,w,[-10,0,10],'k'); hold off
    daspect([1,1,1]); hold off

warning('on', 'MATLAB:contour:ConstantData')
end