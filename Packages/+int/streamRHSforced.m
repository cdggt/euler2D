function bye = streamRHSforced(domain,s,t,loopopts)

D2 = domain.stored.fftD2;
ID2 = domain.stored.fftID2;
Dx = domain.stored.fftDx;
Dy = domain.stored.fftDy;

% T = 7;
% w = 2*pi/T;
x = domain.X-pi;
y = domain.Y-pi;
% k = 2;
% v = domain.nu;
% f = v*sin(x)+sin(k*x).*( (1-2*k*k)*cos(x).*cos(k*y)+4*k*k*k*v*sin(k*y) )./(2*k);
% f = (sin(y).*sin(x)).^2-0.5;
% f = sin(2*y).*sin(2*x);
f = loopopts.F;
% f = cos(x).*cos(y);
% f = util.dealias(domain,f,5);

% x = domain.X;
% y = domain.Y;
% f = sin(2*y);
% f = f-mean(mean(f));
% f = sin(2*x)*2;
f = f*domain.amp;

F = fft2(f);
S = fft2(s);

sx = ifft2( Dx.*S, 'symmetric' );
sy = ifft2( Dy.*S, 'symmetric' );
s2x3y =  ifft2( Dy.*D2.*S, 'symmetric' );
s3x2y =  ifft2( Dx.*D2.*S, 'symmetric' );

bye = fft2(sx.*s2x3y-sy.*s3x2y)-F;
% bye = ifft2(ID2.*bye+domain.nu*D2.*S,'symmetric');
bye = ifft2(ID2.*bye+domain.nu*D2.*S+hv2(S,D2,domain.Nx/4,1/4)+hv2(S,D2,domain.Nx/3,4),'symmetric');

end

% hyperviscosity
function out = hv(S,D2,km,d)
k = sqrt(abs(D2));
% k2 = k(k>km);
% mask = zeros(size(k));
% mask(k>km) = (km-k2).*exp(13./(km-k2));
% out = S.*mask*d;
out = -d*S.*(k.^4)/(km^2);
end
function out = hv2(S,D2,km,d)
k = sqrt(abs(D2));
k2 = k(k>km);
mask = zeros(size(k));
mask(k>km) = (km-k2).*exp(13./(km-k2));
out = S.*mask*d;

end

