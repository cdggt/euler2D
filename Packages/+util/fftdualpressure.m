function [q,term] = fftdualpressure(domain,u,v)
% see int.dualCNRK4qv to see why i pass back the intermediate term. 

    beta = domain.beta;
    Dx = domain.stored.fftDx;
    Dy = domain.stored.fftDy;
    ID2 = domain.stored.fftID2;
    
    ux = u(:,:,1);
    uy = u(:,:,2);
    UX = fft2(ux);
    UY = fft2(uy);
    
    vx = v(:,:,1);
    vy = v(:,:,2);
    VX = fft2(vx);
    VY = fft2(vy);
    
    dxux = real(ifft2(Dx.*UX));
    dyux = real(ifft2(Dy.*UX));
    dxuy = real(ifft2(Dx.*UY));
    dyuy = real(ifft2(Dy.*UY));
    
    dxvx = real(ifft2(Dx.*VX));
    dyvx = real(ifft2(Dy.*VX));
    dxvy = real(ifft2(Dx.*VY));
    dyvy = real(ifft2(Dy.*VY));
    
    udv(:,:,1) = (ux.*dxvx+uy.*dyvx);
    udv(:,:,2) = (ux.*dxvy+uy.*dyvy);
    
    vdu(:,:,1) = (vx.*dxux+vy.*dxuy);
    vdu(:,:,2) = (vx.*dyux+vy.*dyuy);
    
    term = vdu-udv;
    divterm = util.fftdiv(domain,term);
    q = beta*real(ifft2(ID2.*fft2(divterm)));
    
end
