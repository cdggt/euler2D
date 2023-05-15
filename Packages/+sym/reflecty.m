function [state] = reflecty(domain,state)

% Glide reflection symmetry by reflecting over x axis and shifting up half
% a wavelength 

state(:,:,1) =  flip(state(:,:,1),1);
state(:,:,2) = -flip(state(:,:,2),1);


end

