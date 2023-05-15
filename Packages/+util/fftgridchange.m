function [newstate] = fftgridchange(state,sx,sy,dx,dy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the flow field (state) from size (Ni x Ni x 2) to size
% (Nf x Nf x 2). Future iterations will hopefully support non
% square flow fields, but for now only square fields are 
% supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~(exist('dx','var')))
    dx = sy;
    dy = sy;
    sy = sx;
end

normcoeff = (dx/sx)*(dy/sy);
    
if(dx==sx&&dy==sy)
    
    newstate = state;
    
elseif(dx>=sx&&dy>=sy)

    if(mod(sx,2))
        ax = ceil((dx-sx)/2)+1;
        bx = ceil((dx+sx)/2);
    else
        ax = floor((dx-sx)/2)+1;
        bx = floor((dx+sx)/2);
    end
    if(mod(sy,2))
        ay = ceil((dy-sy)/2)+1;
        by = ceil((dy+sy)/2);
    else
        ay = floor((dy-sy)/2)+1;
        by = floor((dy+sy)/2);
    end
    
    for k = 1:size(state,3)
        U = zeros(dy,dx);
        U(ay:by,ax:bx) = fftshift(fft2(state(:,:,k)));
        newstate(:,:,k)  = real(ifft2(ifftshift(U)));  
    end
    
%     UX = zeros(dy,dx);
%     UY = zeros(dy,dx);
%     UX(ay:by,ax:bx) = fftshift(fft2(state(:,:,1)));
%     UY(ay:by,ax:bx) = fftshift(fft2(state(:,:,2)));
% 
%     newstate(:,:,1)  = real(ifft2(ifftshift(UX)));
%     newstate(:,:,2)  = real(ifft2(ifftshift(UY)));  

    newstate = normcoeff.*newstate;
    
elseif(dx>=sx&&dy<sy)

    if(mod(sx,2))
        ax = ceil((dx-sx)/2)+1;
        bx = ceil((dx+sx)/2);
    else
        ax = floor((dx-sx)/2)+1;
        bx = floor((dx+sx)/2);
    end
    if(mod(dy,2))
        ay = ceil((sy-dy)/2)+1;
        by = ceil((sy+dy)/2);
    else
        ay = floor((sy-dy)/2)+1;
        by = floor((sy+dy)/2);
    end
    
    for k = 1:size(state,3)
        Us = fftshift(fft2(state(:,:,k)));
        Ud = zeros(dy,dx);
        Ud(:,ax:bx) = Us(ay:by,:);
        newstate(:,:,k)  = real(ifft2(ifftshift(Ud)));  
    end
    
%     UXs = fftshift(fft2(state(:,:,1)));
%     UYs = fftshift(fft2(state(:,:,2)));
%     
%     UXd = zeros(dy,dx);
%     UYd = zeros(dy,dx);
%     
%     UXd(:,ax:bx) = UXs(ay:by,:);
%     UYd(:,ax:bx) = UYs(ay:by,:);
% 
%     newstate(:,:,1)  = real(ifft2(ifftshift(UXd)));
%     newstate(:,:,2)  = real(ifft2(ifftshift(UYd)));  
    newstate = normcoeff.*newstate;
    
elseif(dx<sx&&dy>=sy)

    if(mod(dx,2))
        ax = ceil((sx-dx)/2)+1;
        bx = ceil((sx+dx)/2);
    else
        ax = floor((sx-dx)/2)+1;
        bx = floor((sx+dx)/2);
    end
    if(mod(sy,2))
        ay = ceil((dy-sy)/2)+1;
        by = ceil((dy+sy)/2);
    else
        ay = floor((dy-sy)/2)+1;
        by = floor((dy+sy)/2);
    end
    
    for k = 1:size(state,3)
        Us = fftshift(fft2(state(:,:,k)));
        Ud = zeros(dy,dx);
        Ud(ay:by,:) = Us(:,ax:bx);
        newstate(:,:,k)  = real(ifft2(ifftshift(Ud)));  
    end
    
%     UXs = fftshift(fft2(state(:,:,1)));
%     UYs = fftshift(fft2(state(:,:,2)));
%     
%     UXd = zeros(dy,dx);
%     UYd = zeros(dy,dx);
%     
%     UXd(ay:by,:) = UXs(:,ax:bx);
%     UYd(ay:by,:) = UYs(:,ax:bx);
% 
%     newstate(:,:,1)  = real(ifft2(ifftshift(UXd)));
%     newstate(:,:,2)  = real(ifft2(ifftshift(UYd)));  
    newstate = normcoeff.*newstate;
    
elseif(dx<sx&&dy<sy)

%     UX = fftshift(fft2(state(:,:,1)));
%     UY = fftshift(fft2(state(:,:,2)));
     
    if(mod(dx,2))
        ax = ceil((sx-dx)/2)+1;
        bx = ceil((sx+dx)/2);
    else
        ax = floor((sx-dx)/2)+1;
        bx = floor((sx+dx)/2);
    end
    if(mod(dy,2))
        ay = ceil((sy-dy)/2)+1;
        by = ceil((sy+dy)/2);
    else
        ay = floor((sy-dy)/2)+1;
        by = floor((sy+dy)/2);
    end
    
    for k = 1:size(state,3)
        U = fftshift(fft2(state(:,:,k)));
        newstate(:,:,k)  = real(ifft2(ifftshift(U(ay:by,ax:bx))));
    end
    
%     newstate(:,:,1)  = real(ifft2(ifftshift(UX(ay:by,ax:bx))));
%     newstate(:,:,2)  = real(ifft2(ifftshift(UY(ay:by,ax:bx))));
    newstate = normcoeff.*newstate;
end

end
