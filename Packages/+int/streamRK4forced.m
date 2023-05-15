function s = streamRK4forced(domain,s,h,t,loopopts)
f = @(s,t) int.streamRHSforced2(domain,s,t,loopopts);

k1 = h*f(s,t);
k2 = h*f(s+0.5*k1,t+h/2);
k3 = h*f(s+0.5*k2,t+h/2);
k4 = h*f(s+k3,t+h);

s = s+( (k1+k4)/2 + k2+k3)/3;

end