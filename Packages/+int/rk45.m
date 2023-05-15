function [p,hh] = rk45(domain,u,t,h,tol,minStep,maxStep)
% h -- initial step size
% tol -- tolerance for determining step size
domain = dom.makeDomain(domain);

T = t(end);
times = zeros(size(t));
hh = zeros(size(t));
hh(1) = h;
t = [t,nan];
f = @(s) int.streamRHSforced2(domain,s);

p = zeros(domain.Ny,domain.Nx,2,length(t)-1);
p(:,:,:,1) = u;
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
    b1 = s2;
    t1 = t2;
    h = h2+1;
    iter = 0;
%     while h/h2 > 1.001 && iter < 2
        iter = iter+1;
        h = h2;
        t2 = t1+h;
        j1 = h*f(b1);
        j2 = h*f(b1+B(2,1)*j1);
        j3 = h*f(b1+B(3,1)*j1+B(3,2)*j2);
        j4 = h*f(b1+B(4,1)*j1+B(4,2)*j2+B(4,3)*j3);
        j5 = h*f(b1+B(5,1)*j1+B(5,2)*j2+B(5,3)*j3+B(5,4)*j4);
        j6 = h*f(b1+B(6,1)*j1+B(6,2)*j2+B(6,3)*j3+B(6,4)*j4+B(6,5)*j5);
        s2 = b1+CH(1)*j1+CH(3)*j3+CH(4)*j4+CH(5)*j5+CH(6)*j6;
        TE = util.fourierMask(domain,CT(1)*j1+CT(3)*j3+CT(4)*j4+CT(5)*j5+CT(6)*j6,mask );
        TE = 512*512*sum(sum(abs( TE )))/numel(b1);
        h2 = tol*0.9*h*( (1/TE)^0.2);
        h2 = max(h2,minStep);
        h2 = min(h2,maxStep);
%     end
%     hh = [hh,h];
%     if isnan(sum(sum(s2)))
%         return
%     end
    if t2 > tt
        hh(tIndex) = h;
        
        
        s = (s2-b1)*(tt-t1)/(t2-t1)+b1;
        p(:,:,:,tIndex) = util.fftunstream(domain,s);
        
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

% u2 = util.fftunstream(domain,s2);






end