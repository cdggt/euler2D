function [] = save_wrapped_omega(title,domain,path,every)

N = size(path,4);
t = linspace(0,domain.Lt,domain.Nt+1);

percent = .2;
bufferx = round(domain.Nx*percent);
buffery = round(domain.Ny*percent);

WriteVideo(title,@(idx) draw(round(every*idx),domain,path,t(round(every*idx)),t(end),bufferx,buffery),floor(N/every));

end

function [] = draw(idx,domain,path,t,T,bufferx,buffery)

omega = util.fftcurl(domain,path(:,:,:,idx));

omega = [omega(:,end-bufferx+1:end) omega omega(:,1:bufferx)];
omega = [omega(end-buffery+1:end,:); omega; omega(1:buffery,:)];
    
colormap parula
contourf(omega, 15,'EdgeColor','None');
hold on
plot([1 domain.Nx+2*bufferx],[buffery buffery],'--k');
plot([1 domain.Nx+2*bufferx],[domain.Ny+buffery domain.Ny+buffery],'--k');
plot([bufferx bufferx],[0 domain.Ny+2*buffery],'--k');
plot([domain.Nx+bufferx domain.Nx+bufferx],[0 domain.Ny+2*buffery],'--k');

[n,m]=size(omega);
xlim([1 n]);
ylim([1 m]);
xticks([]);
yticks([]);
colorbar
title(sprintf('T = %.3f/%.3f s',t,T));

end