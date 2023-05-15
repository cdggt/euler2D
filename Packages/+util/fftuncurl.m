function state = fftuncurl(domain,omega)
    
    Omega = (fft2(omega));

    Psi = -domain.stored.fftID2.*Omega;

    U = domain.stored.fftDy.*Psi;
    V = -domain.stored.fftDx.*Psi;
    state(:,:,1) = real(ifft2((U)));
    state(:,:,2) = real(ifft2((V)));
    
end
