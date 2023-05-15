function [m10,m04,m14] = farazmodes(domain,s,res)
% see 
% A variational approach to probing extreme events in turbulent dynamical systems
% Mohammad Farazmand and Themistoklis P. Sapsis
% for more on why these modes matter in predicting extreme events

if(~exist('res','var'))
    res = 1;
end
n = size(s,4);
indices = 1:res:n;
n = numel(indices);

m10 = zeros(n,1);
m04 = zeros(n,1);
m14 = zeros(n,1);
for j = 1:n
    i = indices(j);
    Omega = fft2(util.fftcurl(domain,s(:,:,:,i)));
    m10(j) = abs(Omega(0+1,1+1));
    m04(j) = abs(Omega(4+1,0+1));
    m14(j) = abs(Omega(4+1,1+1));
end

m10 = 2.*m10./(domain.Nx*domain.Ny);
m04 = 2.*m04./(domain.Nx*domain.Ny);
m14 = 2.*m14./(domain.Nx*domain.Ny);

end