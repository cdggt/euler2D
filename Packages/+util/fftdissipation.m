function D = fftdissipation(domain,s,res)

if(~exist('res','var'))
    res = 1;
end
n = size(s,4);
indices = 1:res:n;
n = numel(indices);

D = zeros(n,1);
Kx = domain.stored.fftDx./(2*pi/domain.Lx);
Ky = domain.stored.fftDy./(2*pi/domain.Ly);
Kn = (Kx.*conj(Kx)+Ky.*conj(Ky));
for j = 1:n
    i = indices(j);
    U = fft2(s(:,:,1,i))./domain.Nx./domain.Ny;
    V = fft2(s(:,:,2,i))./domain.Nx./domain.Ny;
    Di = Kn.*(U.*conj(U)+V.*conj(V));
    D(j) = domain.nu*sum(Di(:));
end

end