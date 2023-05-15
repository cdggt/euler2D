function [state] = dealias(domain,state,kmax)

if(kmax>0&&kmax<inf)
    mask = domain.stored.fftD2;
    mask = (sqrt(abs(mask))<kmax);
    for k = 1:size(state,3)
        state(:,:,k) = ifft2(mask.*fft2(state(:,:,k)));
    end
end

end

