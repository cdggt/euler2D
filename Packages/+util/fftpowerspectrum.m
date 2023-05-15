function [spectrum] = fftpowerspectrum(state)
%GENERATEPOWERSPECTRUM Summary of this function goes here
%   Detailed explanation goes here
[l,w,~] = size(state);
U = fftshift(fft2(state(:,:,1)));
V = fftshift(fft2(state(:,:,2)));
PU = (conj(U).*U)./(l*w);
PV = (conj(V).*V)./(l*w);
spectrum = PU+PV;
warning('remember power spectrum fft is shifted')
end

