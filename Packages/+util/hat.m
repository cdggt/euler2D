function [ uhat ] = hat(domain,u)
    
    norm = util.ffttwonorm(domain,u);
    if(norm>0)
        uhat = u./norm;
    else
        uhat = u;
    end
    
end

