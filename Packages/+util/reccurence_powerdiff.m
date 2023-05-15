function R = reccurence_powerdiff(domain,s,tau,res)

% twomodes = (exist('secondaryy','var')&&exist('secondaryx','var'));
% if(~exist('res','var'))
%     if(~exist('secondaryy','var')&&exist('secondaryx','var'))
%         res = secondaryx;
%     else
%         res = 1;
%     end
% end
if(~exist('res','var'))
    res = 1;
end
n = size(s,4);
ntau = ceil(domain.npu*tau);
indices = 1:res:n;
tauindices = 1:res:ntau;
n = numel(indices);
ntau = numel(tauindices);
R = zeros(numel(tauindices),numel(indices));
for j = 1:n
    i = indices(j);
    s1 = s(:,:,:,i);
    s1 = 2./(domain.Nx*domain.Ny).*abs(fft2(util.fftcurl(domain,s1)));
    for k = 1:ntau
        l = tauindices(k);
        if(i-l>1)
            s2 = s(:,:,:,i-l);
            s2 = 2./(domain.Nx*domain.Ny).*abs(fft2(util.fftcurl(domain,s2)));
            R(k+1,j) = util.ffttwonorm(domain,s1-s2);
        end
    end
end

end