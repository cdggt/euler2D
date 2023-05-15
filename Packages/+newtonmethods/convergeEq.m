function [u,domain] = convergeEq(domain,u,name,tol)


% [xPhase,yPhase] = findPhase(domain,u);
% u = sym.transx(domain,sym.transy(domain,u,yPhase),xPhase);

% define the min_fun and Jacobian from options
n = ceil((domain.Nx/3));
E = util.Energy(u);
    jac = @(domain,b,x,w,h) PO_Jacobian(domain,b,x,w,h,E,n); %...
%          + dot(w,DX(domain,x)) + dot(w,DY(domain,x));
    f = @(x) Min_Fun(domain,x,E,n);

% initialize state for newton to do newton
% domain.kmax = opts.dealias;
w = util.fftcurl(domain,u);
w = util.fftgridchange(w,size(w,1),n);
w = reshape(w, [1,numel(w)])';
w = [w;0;0;0];
b1 = f(w);
fprintf('%.5e',norm(b1)/norm(w));
fprintf('\n');
% newton
iters = 0;
eta = 0.9+rand/11; %opts.eta(3);
q = 0.001;
save_intermediate(domain,w,name,n)
x = linspace(0,2*pi,n+1);
[x,y] = meshgrid(x(1:end-1),x(1:end-1));
r = reshape(b1(1:end-3),[n,n]);
r = util.fftgridchange(r,n,domain.Nx);
r = util.fourierMask(domain,r,@(k) exp(-k/2)./(1+k));
r = util.fftgridchange(r,domain.Nx,n);
r = reshape(r,[n*n,1]);
r = [r;0;0;0];
r = r/norm(r);
r = b1;
while (norm(b1)/norm(w(1:end-3)))>tol && iters<1000
    q = min(1,(q.^0.9)*(1.4+rand/3+1./(1+0.3*iters)));
%     domain.hv = [domain.Nx*31/512,(norm(b1)/norm(w(1:end-1)))*11];
    domain = dom.makeDomain(domain);
    b2 = b1;
    iters = iters+1; 
    % jacobian
    J1 = @(x,w) jac(domain,b1,x,w,1/1234567);
    J = @(q) jac(domain,b1,w,q,1/1234567);

r = b1;
r = J(r);
r = r/norm(r);
    [r,~,N] = newtonmethods.gmres3(J,-b1,r,n*n,eta);
%     [r,~,N] = gmres2(J1,f, w, eta, n*n);
%     toc;
    w1 = w;
    w = w1+q*r;
    b1 = f(w); % update b
        k = 0;
%if (isnan(norm(b1)) || norm(b1)>norm(b2))
        [q,w,b1] = parabolicSearch(f,b1,b2,w1,r,q);
        w = w1+q*r;
        b1 = f(w);
%end
 while (isnan(norm(b1)) || norm(b1)>norm(b2)) && k<77
            k = k+1;
            q = q/(1.1+sqrt(k/10)+rand/11);
            w = w1+q*r;
            b1 = f(w);
            
 end
    
    
    fprintf('%.5e',norm(b1)/norm(w(1:end-3)));
    fprintf(['    ',num2str(N),' iterations','    eta = ', num2str(eta)]);
    fprintf(['    step size ',num2str(q)]);
    fprintf('\n');
    save_intermediate(domain,w,name,n)
    eta = 0.001+rand*0.9+3*0.5/(iters+sqrt(iters)+exp(iters/15));
    eta = min(eta,0.9+rand/11);
end

w = w(1:end-3);
    w = reshape(w,[n,n] );
    w = util.fftgridchange(w,n,domain.Nx);
    u = util.fftuncurl(domain,w);
end


function save_intermediate(domain,s,name,n)
    s = s(1:end-3);
    s = reshape(s,[n,n] );
    s = util.fftgridchange(s,n,domain.Nx);
    u = util.fftuncurl(domain,s);
    save(name,'domain','u');
    w = util.fftcurl(domain,u);
    imagesc(w,[-(max(max(abs(w)))),max(max(abs(w)))]);
    colormap(BWR2())
    daspect([1,1,1])
    set(gca,'YDir','normal')
    pause(0.00001);
end



function b = PO_Jacobian(domain,fx, hi,r,h,E,n)

b = (Min_Fun(domain,hi+h*r,E,n)-fx)/h;

end

function bye = Min_Fun(domain, x,E,n)

Dx = domain.stored.fftDx;
Dy = domain.stored.fftDy;
bye = zeros(size(x));

% reshape, bla bla bla, integrate
w = reshape(x(1:end-3),n,n);
w = util.fftgridchange(w,n,domain.Nx);
u = util.fftuncurl(domain,w);
% [xPhase,yPhase] = findPhase(domain,u);
% w2 = sym.transx(domain,sym.transy(domain,w,yPhase),xPhase);
ux = u(:,:,1);
uy = u(:,:,2);
% Ux = fft2(ux);
% Uy = fft2(uy);
W = fft2(w);



nse = ux.*ifft2(Dx.*W,'symmetric')+uy.*ifft2(Dy.*W,'symmetric');
nse = util.fftgridchange(nse,domain.Nx,n);


% E1 = util.Energy(u1););
bye(1:end-3) = reshape(nse,[n*n,1]);
bye(end) = (E-util.Energy(u))*numel(x)/128;
% bye(end-1) = xPhase;
% bye(end-2) = yPhase;

end

function [q,s,b] = parabolicSearch(f,b1,b2,s,r,q)
qq = [0,0.1+rand/123,0.7+rand*0.2,q];
ff = [norm(b2),norm(f(s+qq(2)*r)),norm(f(s+qq(3)*r)),norm(b1)];
p = polyfit(qq,ff,2);
q = -0.5*p(2)/p(1);
qq = qq(2:end);
ff = ff(2:end);
if q<0 || q>1
    [~,j] = min(ff);
    q = qq(j);
end
b = f(s+q*r);
if norm(b)>norm(b1)
        q = 1;
    b = b1;
    s = s+q*r;
    return
end
s = s+q*r;

end