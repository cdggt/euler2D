function s = streamRK4(domain,s,h)
f = @(s) int.streamRHS(domain,s);

k1 = h*f(s);
k2 = h*f(s+0.5*k1);
k3 = h*f(s+0.5*k2);
k4 = h*f(s+k3);

s = s+( (k1+k4)/2 + k2+k3)/3;

end