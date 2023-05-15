function [X,Y,T,f] = findclosestpass(domain,state,path)

s4 = size(path,4);
F = zeros([domain.Ny,domain.Nx,s4]);

P = fft2(util.hat(domain,util.fftcurl(domain,state)));
for i = 1:size(path,4)
    Q = fft2(util.hat(domain,util.fftcurl(domain,path(:,:,:,i))));
    F(:,:,i) = fft2(conj(Q).*P).*(domain.dx/domain.Nx*domain.dy/domain.Ny);
end

if(s4>1)
% F = squeeze(F(1,:,:));
[percentdistance,idx1d] = max(F(:));
[yidx,xidx,tidx] = ind2sub(size(F),idx1d);
% [yidx,xidx,tidx] = find(F==percentdistance);
X = (min(xidx)-1)*domain.Lx/domain.Nx;
Y = (min(yidx)-1)*domain.Lx/domain.Nx;
T = (min(tidx)-1)*domain.Lt/domain.Nt;
else
% F = squeeze(F(1,:));
[percentdistance,idx1d] = max(F(:));
[yidx,xidx] = ind2sub(size(F),idx1d);
% [yidx,xidx] = find(F==percentdistance);
X = (min(xidx)-1)*domain.Lx/domain.Nx;
Y = (min(yidx)-1)*domain.Lx/domain.Nx;
T = 0; 
end

f = percentdistance;

end