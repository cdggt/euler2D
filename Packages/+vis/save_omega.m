function [] = save_omega(title,domain,path,every)

if(~exist('every','var'))
    every = 1;
end
N = size(path,4);
t = linspace(0,domain.Lt,domain.Nt+1);

WriteVideo(title,@(idx) draw(round(every*idx),domain,path,t(round(every*idx)),t(end)),floor(N/every));

end

function [] = draw(idx,domain,path,t,T)

omega = util.fftcurl(domain,path(:,:,:,idx));
colormap parula
contourf(domain.X,domain.Y,omega, 15,'EdgeColor','None');
colorbar
title(sprintf('T = %.3f/%.3f s',t,T));

end