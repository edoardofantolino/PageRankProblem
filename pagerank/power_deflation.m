function [lambda_2, v_2] = power_deflation(A,lambda_1,v_1,tol,nmax,x0)

%power_deflation(A,lambda_1,v_1,tol,nmax,x0) performs the Wielandt
%deflation to approximate the second most dominant eigenvalue lambda_2 and an
%associated eigenvector v_2 of a n x n matrix A given an approximation 
%lambda_1 to the dominant eigenvalue, an approximation v_1 to the 
%corresponding eigenvector, an absolute error tolerance, a maximum number of 
%iterations and a vector x0 in R^{n-1}

[n,m]= size(A);
if n ~= m
    error ('Only for square matrices'); 
end

if nargin == 3
    x0 = rand (n-1,1);
    nmax = 100;
    tol = 1e-06;
end


i = min(find(abs(v_1) == max(abs(v_1))));


if(i ~= 1)
    B(1:i-1, 1:i-1) = A(1:i-1,1:i-1)-v_1(1:i-1)/v_1(i)*A(i,1:i-1);
    if(i ~= n)
        B(i:n-1,1:i-1) = A(i+1:n,1:i-1)-v_1(i+1:n)/v_1(i)*A(i,1:i-1);
        B(1:i-1,i:n-1) = A(1:i-1,i+1:n)-v_1(1:i-1)/v_1(i)*A(i,i+1:n);
    end
end
if(i ~= n)
    B(i:n-1,i:n-1) = A(i+1:n,i+1:n)-v_1(i+1:n)/v_1(i)*A(i,i+1:n);
end


[lambda_2 ,w_prime, iter]= power_method(B, tol, nmax, x0);

% w_prime = w_prime/norm(w_prime);

if(i ~= 1)
    w(1:i-1,1) = w_prime(1:i-1);
end

w(i) = 0;

if(i ~= n)
    w(i+1:n,1) = w_prime(i:n-1);
end

v_2 = (lambda_2-lambda_1)*w+ (A(i,:)*w)'*(v_1/v_1(i));


end




