function psi = fftstream(domain,state)
    
    Omega = (fft2(util.fftcurl(domain,state)));
    Psi   = -domain.stored.fftID2.*Omega;
    psi   = real(ifft2(Psi));
    
end
