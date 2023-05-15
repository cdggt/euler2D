function [state] = fMask(state,mask)

for k = 1:size(state,3)
    state(:,:,k) = ifft2(mask.*fft2(state(:,:,k)));
end

end