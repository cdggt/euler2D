function vis_fourier_stream(domain,field)

switch size(field,3)
    case 2
        stream = util.fftstream(domain,field);
    case 1
        stream = field;
    otherwise
        warning('field is not velocity or stream function field');
        return;
end

vis.vis_fourier_1dfield(domain,stream,'$\mathcal{F}[\psi]$');

end