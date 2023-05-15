function omega = fftcurl(domain,state)

    Dx = domain.stored.fftDx;
    Dy = domain.stored.fftDy;
    UX = (fft2(state(:,:,1))); 
    UY = (fft2(state(:,:,2)));  
    dyux = real(ifft2((Dy.*UX))); 
    dxuy = real(ifft2((Dx.*UY))); 
    omega = dxuy-dyux; 
end