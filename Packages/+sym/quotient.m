function state = quotient(domain,state,modex,modey,secondaryx,secondaryy)
%GETTEMPLATESTATE
% in this system, there is continuous symmetry about the x axis. To 
% quotient out this symmetry, we use the method of slices ("Reduction of
% SO(2) symmetry for spatially extended dynamical systems" Cvitanovic). 
% This function returns the template point for a state. 

w = util.fftcurl(domain,state);
W = fft2(w);
Fn = W(1+modey,1+modex);
rn = real(Fn);
in = imag(Fn);
mag = abs(Fn/domain.Nx/domain.Ny);
if(mag<1e-8)
    warning('amplitude being quotiented is small')
end

if(rn>0&&in>0)
    phi = atan(abs(in/rn));
elseif(rn<0&&in>0)
    phi = pi-atan(abs(in/rn));
elseif(rn<0&&in<0)
    phi = pi+atan(abs(in/rn));
else
    phi = 2*pi-atan(abs(in/rn));
end

% if(exist('secondaryx','var')&&exist('secondaryy','var'))
%     n = modex^2+modey^2;
%     m = secondaryx^2+secondaryy^2;
%     Fm = W(1+secondaryy,1+secondaryx);
%     rm = real(Fm);
%     im = imag(Fm);
%     phi = norm([angle(Fn),angle(Fm)],2);
% else
%     phi = angle(in/rn);
% end

state = sym.transx(domain,state,-mod(phi,2*pi)+pi);
end

