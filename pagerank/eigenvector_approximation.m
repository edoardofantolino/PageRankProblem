close all
clear
clc

A = [1 2; 3 1];
% teta = pi/6;
% A = [cos(teta) sin(teta); -sin(teta) cos(teta)];
v = [2 0.5]';

[V,L] = eigs(A);

v = v/norm(v)

scatter(V(1,1), V(2,1), 'go', 'filled')
axis square
xlim([0,1])
ylim([0,1])
hold on
grid on
scatter(v(1), v(2), 'ro', 'filled')
pause(7)


for i=1:25
    i
   v = A*v/norm(A*v);
   scatter(v(1), v(2), 'ko')
   pause(1)
end

l = max(A*v)/max(v)

% v2 = A*v/norm(A*v)
% v3 = A^2*v/norm(A^2*v)
% 
% scatter(V(1,1), V(2,1), 'go', 'filled')
% xlim([0,1])
% ylim([0,1])
% hold on
% grid on
% scatter(v(1), v(2), 'ko')
% scatter(v2(1), v2(2), 'ko')
