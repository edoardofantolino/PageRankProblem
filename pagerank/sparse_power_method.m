function[lambda, x, iter] = sparse_power_method(alpha, G, D, e, z, tol, nmax, x0)

A = alpha*G*D;

[n,m] = size(A);
if n~=m
    error('Only for square matrices')
end

if nargin == 5  % number of imput argument, if we select only A
    tol = 1e-8;
    x0 = ones(n,1);
    nmax = 1000;
end

%initialization
x0 = x0/norm(x0);
x = A*x0;
lambda = x0'*x;
err = tol*abs(lambda) + 1;
iter = 0;

for i=1:100
    x = A*x + e*z'*x;
    x = x/norm(x);
    lambda = x'*x;
end

end