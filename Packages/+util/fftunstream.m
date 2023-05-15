function state = fftunstream(domain,psi)
    
    Psi = (fft2(psi));

    U = domain.stored.fftDy.*Psi;
    V = -domain.stored.fftDx.*Psi;
    state(:,:,1) = ifft2((U),'symmetric');
    state(:,:,2) = ifft2((V),'symmetric');
    
end
