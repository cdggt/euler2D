function [state] = gealias(domain,state,kmax,sigma)


mask1 = (sqrt(abs(domain.stored.fftD2))<kmax);
mask2 = mask1<1;
mask = mask1+mask2.*exp(-sigma*sqrt(abs(domain.stored.fftD2)));
for k = 1:size(state,3)
    state(:,:,k) = ifft2(mask.*fft2(state(:,:,k)));
end


end
