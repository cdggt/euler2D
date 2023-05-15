function [state] = transy(domain,state,l)

m = size(state,2);
ky = (-m/2:m/2-1);
T = diag(exp(-1j.*ky*l));

for i = 1:size(state,3)
    U = fftshift(fft2(state(:,:,i)));
    state(:,:,i) = real(ifft2(ifftshift(T*U)));
end

end

