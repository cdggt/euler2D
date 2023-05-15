function [state] = fourierMask(domain,state,mask)

k = sqrt(abs(domain.stored.fftD2));

for j = 1:size(state,3)
    state(:,:,j) = ifft2(mask(k).*fft2(state(:,:,j)),'symmetric');
end

end 