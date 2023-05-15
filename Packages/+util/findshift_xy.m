function [X,Y,f] = findshift_xy(domain,state1,state2)

n1 = util.ffttwonorm(domain,state1);
n2 = util.ffttwonorm(domain,state2);

F = zeros([domain.Ny domain.Nx]);
for i = 1:size(state1,3)
    P = fft2(state1(:,:,i))/n1;
    Q = fft2(state2(:,:,i))/n2;
    F = F + conj(Q).*P;
end
F = real(fft2(F).*(domain.dx/domain.Nx*domain.dy/domain.Ny));

[costheta,idx] = max(F(:));
[yidx,xidx] = ind2sub(size(F),idx);
X = (min(xidx)-1)*domain.Lx/domain.Nx;
Y = (min(yidx)-1)*domain.Ly/domain.Ny;

f = costheta;

end