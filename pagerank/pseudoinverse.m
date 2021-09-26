close all
clear
clc

A = [1 2 3; 2 1 3];

[U,S,V]=svd(A);

Ainv = V*pinv(S)*U'
pinv(A)