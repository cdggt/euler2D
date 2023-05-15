function [state] = dealias2(domain,state,kmax,kmin)

if(kmax>0)
    mask = (sqrt(abs(domain.stored.fftD2))<kmax);
    for k = 1:size(state,3)
        state(:,:,k) = ifft2(mask.*fft2(state(:,:,k)));
    end
end

if(kmin>0)
    mask = (sqrt(abs(domain.stored.fftD2))>kmin);
    for k = 1:size(state,3)
        state(:,:,k) = ifft2(mask.*fft2(state(:,:,k)));
    end
end

end