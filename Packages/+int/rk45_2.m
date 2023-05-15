function [p,u2,hh] = rk45_2(domain,u,t,h,tol,minStep,maxStep,n)
% h -- initial step size
% tol -- tolerance for determining step size
domain = dom.makeDomain(domain);

T = t(end);
times = zeros(size(t));
hh = zeros(size(t));
hh(1) = h;
t = [t,nan];
f = @(s) int.streamRHSforced2(domain,s);

p = zeros(n,n,2,length(t)-1);
p(:,:,:,1) = util.fftgridchange(u,domain.Nx,n);
t2 = 0;
tt = t(2);
tIndex = 2;
s2 = util.fftstream(domain,u);


B = [0, 0, 0 , 0 , 0; ...
    2/9, 0, 0, 0, 0; ...
    1/12, 1/4, 0, 0, 0; ...
    69/128, -243/218, 135/64, 0, 0; ...
    -17/12, 27/4, -27/5, 16/15, 0; ...
    65/432, -5/16, 13/16, 4/27, 5/144];
    
CH = [47/450, 0, 12/25, 32/225, 1/30, 6/25];
CT = [-1/150, 0, 3/100, -16/75, -1/20, 6/25];

hh(1) = h;
h2 = h;
tStart = tic;
mask = @(k) 0.1+k+0.2*k.*k;
while t2 < T
    s1 = s2;
    t1 = t2;
    h = h2+1;
    iter = 0;
%     while h/h2 > 1.001 && iter < 2
        iter = iter+1;
        h = h2;
        t2 = t1+h;
        k1 = h*f(s1);
        k2 = h*f(s1+B(2,1)*k1);
        k3 = h*f(s1+B(3,1)*k1+B(3,2)*k2);
        k4 = h*f(s1+B(4,1)*k1+B(4,2)*k2+B(4,3)*k3);
        k5 = h*f(s1+B(5,1)*k1+B(5,2)*k2+B(5,3)*k3+B(5,4)*k4);
        k6 = h*f(s1+B(6,1)*k1+B(6,2)*k2+B(6,3)*k3+B(6,4)*k4+B(6,5)*k5);
        s2 = s1+CH(1)*k1+CH(3)*k3+CH(4)*k4+CH(5)*k5+CH(6)*k6;
        TE = fMask2(domain,CT(1)*k1+CT(3)*k3+CT(4)*k4+CT(5)*k5+CT(6)*k6,mask );
        TE = 512*512*sum(sum(abs( TE )))/numel(s1);
        h2 = tol*0.9*h*( (1/TE)^0.2);
        h2 = max(h2,minStep);
        h2 = min(h2,maxStep);
%     end
%     hh = [hh,h];
    if isnan(sum(sum(s2)))
        return
    end
    if t2 > tt
        hh(tIndex) = h;
        
        
        s = (s2-s1)*(tt-t1)/(t2-t1)+s1;
        
        
        p(:,:,:,tIndex) = util.fftgridchange(util.fftunstream(domain,s),domain.Nx,n);
        
        if mod(tIndex,round(0.05*(length(t)-2))) == 0
        times(tIndex) = toc(tStart);
        projt = times(tIndex)*(length(t)-2)/(tIndex-1)-times(tIndex);
        hour = floor(projt/3600);
        minn = floor( (projt/60-hour*60) );
        fprintf([num2str(round(100*(tIndex-1)/(length(t)-2))),' percent complete ... estimated ',...
            num2str(hour),' hours and ',num2str(minn),' minutes remaining \n']);
        end
        tIndex = tIndex+1;
        tt = t(tIndex);
    
    end
    
    
end

u2 = util.fftunstream(domain,s2);






end