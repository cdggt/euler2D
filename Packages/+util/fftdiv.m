function divfield = fftdiv(domain, state)

    Dx = domain.stored.fftDx;
    Dy = domain.stored.fftDy;
    U = (fft2(state(:,:,1))); 
    V = (fft2(state(:,:,2)));  
    udx = real(ifft2((Dx.*U))); 
    vdy = real(ifft2((Dy.*V))); 
    divfield = udx+vdy; 
end