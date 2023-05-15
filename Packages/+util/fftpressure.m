function p = fftpressure(domain,state)
    
    beta = domain.beta;
    
    ux = state(:,:,1);
    uy = state(:,:,2);
    UX = fft2(ux);
    UY = fft2(uy);
    Dx = domain.stored.fftDx;
    Dy = domain.stored.fftDy;
    ID2 = domain.stored.fftID2;
    
    dxux = real(ifft2(Dx.*UX));
    dyux = real(ifft2(Dy.*UX));
    dxuy = real(ifft2(Dx.*UY));
    dyuy = real(ifft2(Dy.*UY));
    
    udu(:,:,1) = (ux.*dxux+uy.*dyux);
    udu(:,:,2) = (ux.*dxuy+uy.*dyuy);
    
    term = util.fftdiv(domain,udu);
    p = -beta*real(ifft2(ID2.*fft2(term)));
    
end
