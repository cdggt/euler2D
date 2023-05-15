function [u2,p] = intEuler(domain,u1)
domain = dom.makeDomain(domain);

t = linspace(0,domain.Lt,domain.Nt+1);

h = t(2)-t(1);
s = util.fftstream(domain,u1);
s = util.dealias(domain,s,domain.Nx/3);
if nargout > 1
    p = zeros(domain.Ny,domain.Nx,2,length(t));
    p(:,:,:,1) = u1;
    for j = 2:length(t)
        s = int.streamRK4(domain,s,h);
        u2 = util.fftunstream(domain,s);
        p(:,:,:,j) = u2;
    end
else
for j = 2:length(t)
    s = int.streamRK4(domain,s,h);
end
u2 = util.fftunstream(domain,s);
end

