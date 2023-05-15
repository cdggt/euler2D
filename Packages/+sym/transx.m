function [state] = transx(domain,state,l)

n = size(state,1);
kx = (-n/2:n/2-1);
T = diag(exp(-1j.*kx*l));

for i = 1:size(state,3)
    U = fftshift(fft2(state(:,:,i)));
    state(:,:,i) = real(ifft2(ifftshift(U*T)));
end

end

