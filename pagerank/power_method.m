function[lambda, x, iter] = power_method(A, tol, nmax, x0)

[n,m] = size(A);
if n~=m
    error('Only for square matrices')
end

if nargin == 1  % number of imput argument, if we select only A
    tol = 1e-8;
    x0 = ones(n,1);
    nmax = 1000;
end

%initialization
x0 = x0/norm(x0);
pro = A*x0;
lambda = x0'*pro;
err = tol*abs(lambda) + 1;
iter = 0;

while (err>tol*abs(lambda)) && (iter <= nmax)
    x = pro;
    x = x/norm(x);
    pro = A*x;
    lambdanew = x'*pro;
    err = abs(lambdanew-lambda);
    lambda = lambdanew;
    iter = iter + 1;
end

end