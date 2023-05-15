function [ norm ] = ffttwonorm(domain,u,dx,dy)
    
    if(~exist('dx','var'))&&(isfield(domain,'dx'))
        dx = domain.dx;
    end
    if(~exist('dy','var'))&&(isfield(domain,'dy'))
        dy = domain.dy;
    end
    
    norm = 0;
    [n,m,d] = size(u);
    for i=1:d
        field = u(:,:,i);
        Field = fft2(field);
        magnitudes = Field.*conj(Field);
        norm = norm + sum(magnitudes(:));
    end
    norm = real(sqrt(dx*dy*norm/n/m));
    
end

