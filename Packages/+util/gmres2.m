function [x, error, total_iters] = gmres2(J,f, xc, tol, maxiter, r,f0)
% GMRES linear equation solver for use in Newton-GMRES solver
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any 
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%            params(3) (Optional) = reorthogonalization method
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarted GMRES
%              will be used as the linear solver.
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires givapp.m, dirder.m 

%
% initialization
%
% f0 = f(xc);

reorth = 1;
% if length(params) == 3
%     reorth = params(3);
% end
%
% The right side of the linear equation for the step is -f0. 
%
b = -f0;
n = length(b);
%
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
%
x = zeros(n,1);
% r = b;
% r = r+(rand(size(r))-0.5)*norm(r)/numel(r);
% if nargin == 5
%     x = xinit;
%     r = -dirder(xc, x, f, f0)-f0;
% end
%
%
h = zeros(maxiter);
v = zeros(n,maxiter);
c = zeros(maxiter+1,1);
s = zeros(maxiter+1,1);
rho = norm(r);
g = rho*eye(maxiter+1,1);
tol = tol*norm(b);
error = [];
%
% Test for termination on entry.
%
error = [error,rho];
total_iters = 0;
%
%
v(:,1) = r/rho;
beta = rho;
k = 0;
%
% GMRES iteration
%
while (((rho > tol) && (k < maxiter)) ) || k<1
    k = k+1;
%
%   Call directional derivative function.
%
    v(:,k+1) = J(xc, v(:,k));
    normav = norm(v(:,k+1));
%
%   Modified Gram-Schmidt
%
    for j = 1:k
        h(j,k) = v(:,j)'*v(:,k+1);
        v(:,k+1) = v(:,k+1)-h(j,k)*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
    normav2 = h(k+1,k);
%
%   Reorthogonalize?
%
if  (reorth == 1 && normav + .001*normav2 == normav) || reorth ==  3
    for j = 1:k
        hr = v(:,j)'*v(:,k+1);
	h(j,k) = h(j,k)+hr;
        v(:,k+1) = v(:,k+1)-hr*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
end
%
%   Watch out for happy breakdown.
%
    if(h(k+1,k) ~= 0)
    v(:,k+1) = v(:,k+1)/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation.
%
    if k > 1
        h(1:k,k) = givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
%
%   Don't divide by zero if solution has  been found.
%
    nu = norm(h(k:k+1,k));
    if nu ~= 0
%        c(k) = h(k,k)/nu;
        c(k) = conj(h(k,k)/nu);
        s(k) = -h(k+1,k)/nu;
        h(k,k) = c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k) = 0;
        g(k:k+1) = givapp(c(k),s(k),g(k:k+1),1);
    end
%
%   Update the residual norm.
%
    rho = abs(g(k+1));
    error = [error,rho];
%
%   end of the main while loop
%
end
%
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
%
y = h(1:k,1:k)\g(1:k);
total_iters = k;
x = x + v(1:n,1:k)*y;
end


% function z = dirder(x,w,f,f0)
% % Finite difference directional derivative
% % Approximate f'(x) w
% % 
% % C. T. Kelley, April 1, 2003
% %
% % This code comes with no guarantee or warranty of any kind.
% %
% % function z = dirder(x,w,f,f0)
% %
% % inputs:
% %           x, w = point and direction
% %           f = function
% %           f0 = f(x), in nonlinear iterations
% %                f(x) has usually been computed
% %                before the call to dirder
% 
% %
% % Use a hardwired difference increment.
% %
% epsnew = 1.d-7;
% %
% n = length(x);
% %
% % scale the step
% %
% if norm(w) == 0
%     z = zeros(n,1);
% return
% end
% %
% % Now scale the difference increment.
% %
% xs=(x'*w)/norm(w);
% if xs ~= 0.d0
%      epsnew=epsnew*max(abs(xs),1.d0)*sign(xs);
% end
% epsnew=epsnew/norm(w);
% %
% % del and f1 could share the same space if storage
% % is more important than clarity.
% %
% del = x+epsnew*w;
% f1 = feval(f,del);
% z = (f1 - f0)/epsnew;
% end

function vrot = givapp(c,s,vin,k)
%  Apply a sequence of k Givens rotations, used within gmres codes.
% 
%  C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
%  function vrot = givapp(c, s, vin, k)
%
vrot = vin;
for i = 1:k
    w1 = c(i)*vrot(i)-s(i)*vrot(i+1);
%
%   Here's a modest change that makes the code work in complex
%   arithmetic. Thanks to Howard Elman for this.
%
%    w2 = s(i)*vrot(i)+c(i)*vrot(i+1);
    w2 = s(i)*vrot(i)+conj(c(i))*vrot(i+1);
    vrot(i:i+1) = [w1,w2];
end
end