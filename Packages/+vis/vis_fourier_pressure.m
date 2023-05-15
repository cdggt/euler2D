function vis_fourier_pressure(domain,field)

switch size(field,3)
    case 2
        pressure = util.fftpressure(domain,field);
    case 1
        pressure = field;
    otherwise
        warning('field is not velocity or pressure field');
        return;
end

vis.vis_fourier_1dfield(domain,pressure,'$\mathcal{F}[p]$');

end