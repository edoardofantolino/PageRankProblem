close all
clear
clc
 
A = [2     3     1
     0     0     2
     0     1     2];

% A = [2     3     1
%      3     3     2
%      1     2     2];
 
% [Q1,R1] = qr(A);
% A2 = R1*Q1;
% [Q2,R2] = qr(A2);
% A3 = R2*Q2;
% [Q3,R3] = qr(A3);
% A4 = R3*Q3;
[Q,R] = qr(A)
An = A;
for i=1:15
   [Qn,Rn] = qr(An);
   An = Rn*Qn;
   Q=Qn*Q;
end

[V,L] = eigs(A);
diag(L)
sort(diag(Q'*A*Q), 'descend')