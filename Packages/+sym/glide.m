function [state] = glide(domain,state,n)

% Glide reflection symmetry by reflecting over x axis and shifting up half
% a wavelength 

wavelength = domain.Ly/domain.Ky;
state(:,:,1) = -flip(state(:,:,1),2);
state(:,:,2) =  flip(state(:,:,2),2);
state = sym.transy(domain,state,wavelength/2);

if(exist('n','var')&&n>1)
    state = sym.glide(domain,state,n-1);
end

end

