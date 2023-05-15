function [ norm ] = fftinnerproduct(domain,u,v,dx,dy)
    
    if(~exist('dx','var'))&&(isfield(domain,'dx'))
        dx = domain.dx;
    end
    if(~exist('dy','var'))&&(isfield(domain,'dy'))
        dy = domain.dy;
    end
    
    
    if(size(u)~=size(v))
        error('vectors must be of same dimension');
    end
    [n,m,d] = size(u);
    [~,~,dother] = size(u);
    if(d~=dother)
        warning('you are taking the inner project of states with different dimension.')
    end
    
    norm = 0;
    for i =1:d
        magnitudes = conj(fft2(u(:,:,i))).*fft2(v(:,:,i));
        norm = norm + sum(magnitudes(:));
    end
    norm = real(dx*dy*norm/n/m);
    
end

