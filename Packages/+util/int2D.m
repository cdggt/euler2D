function val = int2D(domain,array,dx,dy)
    
    if(~exist('dx','var'))&&(isfield(domain,'dx'))
        dx = domain.dx;
    end
    if(~exist('dy','var'))&&(isfield(domain,'dy'))
        dy = domain.dy;
    end
    val = dx*dy*trapz(trapz(array));

end