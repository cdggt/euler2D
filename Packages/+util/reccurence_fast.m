function R = reccurence_fast(domain,p,ntau,ntheta,res)

if(~exist('res','var'))
    res = 1;
end
ntau = round(ntau);
n = size(p,4);
% ntau = round(size(p,4));
indices = 1:res:n;
tauindices = 1:res:ntau;
n = numel(indices);
ntau = numel(tauindices);
theta = 0;
% theta = theta(1:end-1);

R = zeros(numel(tauindices),numel(indices));
for j = n:-1:1
    i = indices(j);
    for k = 1:ntau
        l = tauindices(k);
        if(i-l>1)
            s1 = p(:,:,:,i)/sum(sum(sum(abs(p(:,:,:,i)))));
            s2 = p(:,:,:,i-l)/sum(sum(sum(abs(p(:,:,:,i-1)))));
%             fmin = inf;
            f = sum(sum(sum( abs(s1-s2))));
            R(k+1,j) = f;
        end
    end
%     fprintf([num2str(j),' out of ', num2str(n), '\n'])
end

end