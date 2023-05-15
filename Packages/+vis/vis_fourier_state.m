function vis_fourier_state(domain,field)

switch size(field,3)
    case 2
        state = field;
    case 1
        state = util.fftuncurl(domain,field);
    otherwise
        warning('field is not velocity or vorticity field');
        return;
end

vis.vis_fourier_2dfield(domain,state,'$\mathcal{F}[u_x]$','$\mathcal{F}[u_y]$');

end