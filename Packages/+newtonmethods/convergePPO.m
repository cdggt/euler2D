function [u,domain,iters] = convergePPO(domain,u,name,tol,maxiter)




% define the min_fun and Jacobian from options
    jac = @(domain,b,x,w,h) PO_Jacobian(domain,b,x,w,h); %...
%          + dot(w,DX(domain,x)) + dot(w,DY(domain,x));
    f = @(x) Min_Fun_PO(domain,x);

% initialize state for newton to do newton
% domain.kmax = opts.dealias;
s = util.fftcurl(domain,u);
s = reshape(s, [1,numel(s)])';
b1 = f(s);
fprintf('%.5e',norm(b1)/norm(s));
fprintf('\n');
% newton
iters = 0;
eta = 0.2+rand/11; %opts.eta(3);
q = 0.3;
r = b1;
while (norm(b1)/norm(s(1:end-1)))>tol && iters<maxiter
    q = min(1,0.1+(q.^0.9)*(1.1+rand/3));
    b2 = b1;
    iters = iters+1; 
    % jacobian
%     J = @(x,w) jac(domain,b1,x,w,1/123456789);
    J = @(q) jac(domain,b1,s,q,1/1234567);
    r = s;
    r = r/norm(r);
    eta = 0.45+rand*0.54;
    [r,~,n] = newtonmethods.gmres3(J,-b1,r,10,eta);
%     [r,~,n] = gmres2(J,f,s,eta1,777);
%     r = reshape(r,[domain.Ny,domain.Nx]);
%     r = util.dealias(domain,r,16);
%     r = reshape(r,size(w));
ss = s;
    s = s+q*r;
    b1 = f(s); % update b
    bu = util.fftuncurl(domain,reshape(b1,[domain.Nx,domain.Ny]));
    u = util.fftuncurl(domain,reshape(s,[domain.Nx,domain.Ny]));
    zeta_u = sqrt(util.Energy(bu)/util.Energy(u));
    
    k = 0;
if (isnan(norm(b1)) || norm(b1)>norm(b2))
        [q,s,b1] = parabolicSearch(f,b1,b2,ss,r,q);
%         s = ss+q*r;
%         b1 = f(s);
end



 while (isnan(norm(b1)) || norm(b1)>norm(b2)) && k<18
            k = k+1;
            q = q/(1+rand/2+sqrt(k)/4);
            s = ss+q*r;
            dammp23 = 0.2*norm(b1)/norm(s);
            ss = reshape(ss,[domain.Ny,domain.Nx]);
            s = reshape(s,[domain.Ny,domain.Nx]);
%             for j = 1:size(s,3)
%                 ss(:,:,j) = util.gealias(domain,ss(:,:,j),1,dammp23);
%                 s(:,:,j) = util.gealias(domain,s(:,:,j),1,dammp23);
%             end
            s = reshape(s, [numel(s),1]);
            ss = reshape(ss, [numel(s),1]);
            b1 = f(s);
            
end
    
    
    fprintf('%.5e',norm(b1)/norm(s));
    fprintf('    ');
    fprintf('%.5e',zeta_u);
    fprintf(['    ',num2str(n),' iterations','    eta = ', num2str(eta)]);
    fprintf(['    step size ',num2str(q)]);
    fprintf('\n');
    
     % update eta
     
%     if norm(b1) < norm(b2)
%         eta1 = (eta1*10+(1-norm(b1)/norm(b2))^(1/50))/11;
%         eta1=eta1*0;
%     else
%         eta1 = eta1*(1-sqrt(eps))-sqrt(eps);
%     end
%             eta1 = min(eta1,opts.eta(1));
%         eta1 = max(eta1,opts.eta(2));
    save_intermediate(domain,s,name)
    eta = 0.001+rand*0.9+3*0.5/(iters+sqrt(iters)+exp(iters/15));
    eta = min(eta,0.9+rand/11);
end

save_intermediate(domain,s,name)
end

function [q,s,b] = parabolicSearch(f,b1,b2,s,r,q)
if q > 0.9999
    qq = [0,0.1+rand/123,0.7+rand*0.2,q];
else
    qq = [0,0.45+rand/10,1,q];
end
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


function save_intermediate(domain,w,name)
    w = reshape(w,[domain.Ny,domain.Nx]);
    u = util.fftuncurl(domain,w);
        imagesc(w,[-max(max(abs(w))),max(max(abs(w)))]); hold on
    colormap(BWR2())
    colorbar
    set(gca,'YDir','normal')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
        set(gca,'ytick',[])
    set(gca,'yticklabel',[])
         hold off
    daspect([1,1,1]);
    pause(0.00000015);
    save(name,'domain','u');
end



function b = PO_Jacobian(domain,fx, hi,r,h)

b = (Min_Fun_PO(domain,hi+h*r)-fx)/h;

end

function bye = Min_Fun_PO(domain, x)
% initialize output


% reshape, bla bla bla, integrate
w1 = reshape(x,[domain.Ny,domain.Nx,1]);
u1 = util.fftuncurl(domain,w1);
t = [0,domain.Lt];
% [p2,~] = int.rk45(domain,u1,t,1/2000,1.7,1/12345,1/50);
% u2 = p2(:,:,:,end);
u2 = int.intEuler(domain,u1);
w2 = zeros(size(w1)+1);
w2(1:end-1,1:end-1) = util.fftcurl(domain,u2);
w2(end,1:end) = w2(1,1:end);
w2(1:end,end) = w2(1:end,1);
w2 = rot90(w2,-1);
w2 = w2(1:end-1,1:end-1);
u2 = util.fftuncurl(domain,w2);
[xPhase,yPhase] = util.findPhase(domain,u2);
u2 = sym.transx(domain,sym.transy(domain,u2,yPhase),xPhase);
w2 = util.fftcurl(domain,u2);
y = reshape(w2,size(x));

E1 = util.Energy(u1);
E2 = util.Energy(u2);
bye = (x-y*sqrt(E1/E2));
% bye = x-y;
% bye = reshape(fMask2(domain,w1/E1-w2/E2,@(k) 1./(1+77*n*k)),size(x));
 



end