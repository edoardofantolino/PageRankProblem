function[lambda, x, iter] = no_mm_power_method(alpha,G,D,e,z,c_j, tol, nmax, x0)
A = alpha*G*D;

[n,m] = size(A);
if n~=m
    error('Only for square matrices')
end

if nargin == 6  % number of imput argument, if we select only A
    tol = 1e-8;
    x0 = ones(n,1);
    nmax = 1000;
end

%initialization
x0 = x0/norm(x0);

v = x0;
Mnv = zeros(length(G),1);
for j=1:length(G)
    if c_j(j)~=0
        Mnv = Mnv + alpha*G(:,j)*D(j,j)*v(j);
    end
end   

x = Mnv;
lambda = x0'*x;
err = tol*abs(lambda) + 1;
iter = 0;

for i=1:100
    
    v = x;
    Mnv = zeros(length(G),1);
    for j=1:length(G)
        if c_j(j)~=0
            Mnv = Mnv + alpha*G(:,j)*D(j,j)*v(j);
        end
    end   
    
    x = Mnv + e*z'*x;
    
    x = x/norm(x);
    lambda = x'*x;
end

end