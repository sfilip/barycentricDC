% Specify the function to approximate
f = @(x) (0.5*erf(x./sqrt(0.0002))+1.5).*exp(-x);
% Specify the approximation set
t = linspace(-1,1,20000)';
ft = f(t);

% Construct a type (16,16) best approximation
% to this function
% Our tests used the Mosek optimizer
% if available, it can be set by calling
% cvx_solver Mosek
[rh,xk,err] = barycentricDC(ft,t,16,16,'display','log');
