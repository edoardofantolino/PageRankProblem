close all
clear
clc

data = dlmread('toy_example_matlab.txt');
G = sparse(data(:, 1), data(:, 2), 1);

% (*1)
% figure(1)
% spy(G)

% (*3)
out_degree = sum(G);
in_degree = sum(G');

r_j = in_degree;
c_j = out_degree;



% (*5)
% M = G(:,1)/c_j(1);
% for k=2:length(G)
%     if c_j(k)~=0 
%         M = [M G(:,k)/c_j(k)];
%     else 
%         M = [M zeros(length(G), 1)];
%     end
% end
% 
% figure(2)
% spy(M)

% (*7)
Mt = G(:,1)/c_j(1);
for k=2:length(G)
    if c_j(k)~=0 
        Mt = [Mt G(:,k)/c_j(k)];
    else 
        k;
        Mt = [Mt (1/length(G))*ones(length(G),1)];
    end
end

% figure(3)
% spy(Mt)
% 
% lambdas = eigs(Mt);
% 
% power_method(Mt);

% (*12)
% t = (1/length(G)) * ones(length(G),1);
t = zeros(length(G),1);
t(1) = 1;
alpha = 0.85;
delta = (1-alpha)/length(G);
A = alpha*G(:,1)/c_j(1)+delta;
for k=2:length(G)
    if c_j(k)~=0 
        A = [A alpha*G(:,k)/c_j(k)+delta];
    else 
        k;
        A = [A (1/length(G))*ones(length(G),1)];
    end
end

% % (**14)
% 
% [V,L] = eigs(A);
% 
% [lambda, x, iter] = power_method(A);
% 
% x = x/norm(x);
% xMATLAB = V(:,1)/norm(V(:,1));
% 
% rel_err = norm(xMATLAB - x) / norm(xMATLAB);
% 
% (**15)

% A = alpha*G*D + e*z'
alpha = 0.85;
delta = (1-alpha)/length(G);
e = ones(length(G),1);
z = (1/length(G))*ones(length(G),1);

d_j = zeros(length(G),1);
for k=1:length(G)
    if c_j(k) ~= 0
        d_j(k) = 1/c_j(k);
        z(k) = delta;
    end
end

D = spdiags(d_j,0,length(G),length(G));

MATLAB = eigs(A,1);
POWER_METHOD = power_method(A);
[lambda, x, iter] = sparse_power_method(alpha, G, D, e, z);


tic;
[lambda, x, iter] = power_method(A);
toc;
tic;
[lambda, x, iter] = sparse_power_method(alpha, G, D, e, z);
toc;
% % tic;
% % [lambda, x, iter] = bad_sparse_power_method(alpha, G, D, e, z);
% % toc;
% 
% % (**17)
% 
% % [V,L] = eigs(Mt);
% % V(:,1);
% % L(1);
% % 
% % [lambda, x, iter] = power_method(Mt);
% % lambda;
% % x;
% % 
% % [lambda_2, x_2] = power_deflation(Mt,lambda,x)
% % V(:,2)
% % max(L)
% 
% % (**19)
% 
% podium = maxk(x,3)
% 
% for k=1:3
%    find(x==podium(k))
% end