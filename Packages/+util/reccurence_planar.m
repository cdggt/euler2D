function [R,L] = reccurence_planar(domain,s,N)

M = size(s,4);
R = zeros(N,M);
if(nargout>1)
    L = zeros(N,M);
end

for i = 1:M
    s1 = s(:,:,:,i);
    w1 = util.hat(domain,util.fftcurl(domain,s1));
    W1 = fft2(w1);
    for j = 0:N-1
        if(i-j>0)
            s2 = s(:,:,:,i-j);
            w2 = util.hat(domain,util.fftcurl(domain,s2));
            W2 = fft2(w2);
            
            fl = fft2(conj(W2).*W1).*(domain.dx/domain.Nx*domain.dy/domain.Ny);
            fl = fl(1,:);
            [fmax,idx] = max(fl(:));
            lmin = (idx-1)*domain.Lx/domain.Nx;
            
            if(nargout>1)
                L(j+1,i) = lmin;
            end
            R(j+1,i) = 1-fmax;
        end
    end
end
end

% function R = reccurence_spiral(domain,s,tau,res)
% 
% if(~exist('res','var'))
%     res = 1;
% end
% n = size(s,4);
% ntau = ceil(domain.npu*tau);
% indices = 1:res:n;
% tauindices = 1:res:ntau;
% n = numel(indices);
% ntau = numel(tauindices);
% R = zeros(numel(tauindices),numel(indices));
% L = linspace(0,2*pi,domain.Nx+1);
% L = L(1:end-1);
% for j = 1:n
%     i = indices(j);
%     s1 = s(:,:,:,i);
%     w1 = util.fftcurl(domain,s1);
%     for k = 1:ntau
%         l = tauindices(k);
%         if(i-l>1)
%             s2 = s(:,:,:,i-l);
%             w2 = spiral(domain,s2,L);
% %             vis.vis_1dfield(domain,w2);
% %             vis.vis_omega(domain,s2);
%             fx = sum((w1-w2).^2,2);
%             [~,idx] = min(fx);
%             Lmin = L(idx);
%             s2 = sym.transx(domain,s2,Lmin);
%             R(k+1,j) = util.ffttwonorm(domain,s1-s2);
%         end
%     end
% end
% end
% 
% function s = spiral(domain,s,l)
% 
% w = util.fftcurl(domain,s);
% for i = 1:numel(l)
%     w(i,:) = sym.transx(domain,w(i,:),l(i));
% end
% % W = fftshift(fft2(w));
% % n = size(s,1);
% % kx = (-n/2:n/2-1);
% % T = diag(exp(-1j.*(l.^2.*kx)));
% % s = ifft(ifftshift(W.*T),'symmetric');
% s = w;
% 
% end




