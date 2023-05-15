function [newstate] = fftsquaregridchange(state,Ni,Nf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the flow field (state) from size (Ni x Ni x 2) to size
% (Nf x Nf x 2). Future iterations will hopefully support non
% square flow fields, but for now only square fields are 
% supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normcoeff = (Nf/Ni)^2;
isdodd = mod(Nf,2);
issodd = mod(Ni,2);
    
if(Nf==Ni)
    
    newstate = state;
    
elseif(Nf>Ni)

    if(issodd)
        a = ceil((Nf-Ni)/2)+1;
        b = ceil((Nf+Ni)/2);
    else
        a = floor((Nf-Ni)/2)+1;
        b = floor((Nf+Ni)/2);
    end
    
    UX = zeros(Nf,Nf);
    UY = zeros(Nf,Nf);
    UX(a:b,a:b) = fftshift(fft2(state(:,:,1)));
    UY(a:b,a:b) = fftshift(fft2(state(:,:,2)));

    newstate(:,:,1)  = real(ifft2(ifftshift(UX)));
    newstate(:,:,2)  = real(ifft2(ifftshift(UY)));  
    newstate = normcoeff.*newstate;
    
else

    UX = fftshift(fft2(state(:,:,1)));
    UY = fftshift(fft2(state(:,:,2)));
    
    if(isdodd)
        a = ceil((Ni-Nf)/2)+1;
        b = ceil((Ni+Nf)/2);
    else
        a = floor((Ni-Nf)/2)+1;
        b = floor((Ni+Nf)/2);
    end
    
    newstate(:,:,1)  = real(ifft2(ifftshift(UX(a:b,a:b))));
    newstate(:,:,2)  = real(ifft2(ifftshift(UY(a:b,a:b))));
    newstate = normcoeff.*newstate;
end

end
