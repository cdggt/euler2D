function vis_fourier_omega(domain,field)

switch size(field,3)
    case 2
        omega = util.fftcurl(domain,field);
    case 1
        omega = field;
    otherwise
        warning('field is not velocity or vorticity field');
        return;
end

vis.vis_fourier_1dfield(domain,omega,'$\mathcal{F}[\omega]$');

end