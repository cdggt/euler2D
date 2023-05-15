function u = streamFunctionSmooth(domain,blurradius,u)

s = util.fftstream(domain,u);
w = util.fftcurl(domain,u);
n = size(u,1);
wnew = w;
x=domain.X(1,:);
y=domain.Y(:,1);

for i = 1:n
    for j = 1:n
        f = exp(-(s-s(i,j)).^2/2/blurradius^2);
        wnew(i,j) = trapz(y,trapz(x,f.*w,2),1)/trapz(y,trapz(x,f,2),1);
    end
end

u = util.fftuncurl(domain,wnew);

end