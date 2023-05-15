function S = statespacespeed(domain,s,res)

if(~exist('res','var'))
    res = 1;
end
n = size(s,4);
indices = 1:res:n;
n = numel(indices);

S = zeros(n,1);
for j = 1:n
    i = indices(j);
    S(j) = util.ffttwonorm(domain,int.taninf(domain,s(:,:,:,i)))/util.ffttwonorm(domain,s(:,:,:,i));
end

end