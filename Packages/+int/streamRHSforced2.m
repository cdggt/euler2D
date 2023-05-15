function bye = streamRHSforced2(domain,s)

D2 = domain.stored.fftD2;
ID2 = domain.stored.fftID2;
Dx = domain.stored.fftDx;
Dy = domain.stored.fftDy;

% T = 7;
% w = 2*pi/T;
x = domain.X-pi;
y = domain.Y-pi;



F = domain.stored.fftcurlf;
S = fft2(s);

sx = ifft2( Dx.*S, 'symmetric' );
sy = ifft2( Dy.*S, 'symmetric' );
s2x3y =  ifft2( Dy.*D2.*S, 'symmetric' );
s3x2y =  ifft2( Dx.*D2.*S, 'symmetric' );

bye = fft2(sx.*s2x3y-sy.*s3x2y)-F;
% bye = ifft2(ID2.*bye+domain.nu*D2.*S,'symmetric');
k = sqrt(abs(D2));
k(k>domain.Nx*2/6) = 0;
k(k>0) = 1;
bye = ifft2(k.*(ID2.*bye+hv2(S,D2,domain.hv(1),domain.hv(2))+domain.nu*D2.*S),'symmetric')-domain.alpha*s;

end

% hyperviscosity

function out = hv2(S,D2,km,d)
k = sqrt(abs(D2));
k2 = k(k>km);
mask = zeros(size(k));
mask(k>km) = (km-k2).*exp(13./(km-k2));
out = S.*mask*d;

end