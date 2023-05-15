function X = find_phase_x(domain,u1,u2)

% This does the same stuff as find_shift_x but then also uses newton
% iterations to converge on the best shift. The resulution of find_shift_x
% is L_x/N_x. The resolution here is machine precision. 

% josh step
F = zeros([domain.Ny,domain.Nx,1]);

P = fft2(util.hat(domain,util.fftcurl(domain,u1)));
for i = 1:size(u2,4)
    Q = fft2(util.hat(domain,util.fftcurl(domain,u2(:,:,:,i))));
    F(:,:,i) = fft2(conj(Q).*P).*(domain.dx/domain.Nx*domain.dy/domain.Ny);
end

F = F(1,:);
[~,xidx] = max(F(:));
% [~,xidx] = ind2sub(size(F),idx1d);
X = (min(xidx)-1)*domain.Lx/domain.Nx;

% dz step
h = 10^-4.1;
f = @(x) norm( reshape(u2-sym.transx(domain,u1,x),[numel(u1),1]) );
dfdx = @(p) (f(p+h)-f(p-h))/10000;
X = bi_secant(dfdx,X-0.1,X+0.1,10^-12,15,15);
 

end

function [x] = bi_secant(f,x1,x2,tol,iter1,iter2)
x = (x1+x2)/2;
j = 0;
f1 = f(x1);
fm = tol+1;
% bijection
while abs(fm) > tol && j<iter1
    fm = f(x);
    j = j+1;
    if sign(f1) == sign(fm)
        x1 = x;
        f1 = fm;
    else
        x2 = x;
    end
    x = (x1+x2)/2;
end

j=0;
while abs(fm) > tol && j<iter2
    fm = f(x);
    j = j+1;
    df = (f1-fm)/(x1-x);
    f1 = fm;
    x1 = x;
    x = x-fm/df;
end
end