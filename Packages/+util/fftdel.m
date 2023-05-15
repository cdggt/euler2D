function delfield = fftdel(domain,field)

    Dx = domain.stored.fftDx;
    Dy = domain.stored.fftDy;
    F = fft2(field); 
    delfield(:,:,1) = real(ifft2((Dx.*F))); 
    delfield(:,:,2) = real(ifft2((Dy.*F))); 
    
end