function R = reccurence_avgx(domain,s,tau,res)

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
    s1 = avgx(domain,s1);
    for k = 1:ntau
        l = tauindices(k);
        if(i-l>1)
            s2 = s(:,:,:,i-l);
            s2 = avgx(domain,s2);
            R(k+1,j) = util.ffttwonorm(domain,s1-s2);
        end
    end
end

end

function s = avgx(domain,s)

    w = util.fftcurl(domain,s);
    W = fft2(w);
    W(:,2:end) = 0;
    w = ifft2(W,'symmetric');
    s = util.fftuncurl(domain,w);
    
end





