function [s,tau,shift,f] = getminreccurence(domain,path,reccurencematrix,percentcutoff,resolution)

if(~exist('resolution','var'))
    resolution = 1;
end

% find min value in reccurence
ncutoff = round(percentcutoff*size(reccurencematrix,1));
subset = reccurencematrix;
subset(1:ncutoff,:)=inf;
subset(subset==0)=inf;
f = min(subset(:));
[row,col] = find(subset==f);

tau = (row-1)*resolution/domain.npu;
% t = (col-row+1)*resolution/domain.npu;
s = path(:,:,:,col-row+1);
[fcheck,shift] = test_s(domain,s,tau);

end

function [f,L] = test_s(domain,s1,tau)


s2 = int.tangent(domain.append('lt',tau),s1);
s2 = s2(:,:,:,end);

w1 = util.hat(domain,util.fftcurl(domain,s1));
W1 = fft2(w1);
w2 = util.hat(domain,util.fftcurl(domain,s2));
W2 = fft2(w2);
            
fl = fft2(conj(W2).*W1).*(domain.dx/domain.Nx*domain.dy/domain.Ny);
fl = fl(1,:)-1;
[fmin,mindx] = min(abs(fl));
L = (mindx-1)*domain.Lx/domain.Nx;
   
f = fmin;
% f = util.ffttwonorm(domain,s1-sym.transx(domain,s2,-L))/util.ffttwonorm(domain,s1);

end
