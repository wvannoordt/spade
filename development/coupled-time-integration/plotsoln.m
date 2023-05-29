clear
clc
close all

mat = [0.01 -0.05; 0.04 0.03];
[v, lam] = eig(mat);

dt = 1e-2;
n = 20000;
t = linspace(0, (n-1)*dt, n)';
x0 = [1.0;0.5];

sol = zeros(n, 2);

for i=1:n

  sol(i,:) = (e^(t(i)*mat))*x0;

end

data = dlmread('soln.dat');
tnum = data(:,1);
tana = t;
num  = data(:, 2:3);

figure('position', [0 0 1800 800])
subplot(1, 2, 1)
plot(sol(:,1), sol(:,2), 'linewidth', 4);
hold on
plot(num(:,1), num(:,2), 'linewidth', 3);
subplot(1, 2, 2)
plot(tana, sol(:, 1), 'linewidth', 4, 'color', [0 0 0.5]);
hold on
plot(tana, sol(:, 2), 'linewidth', 4, 'color', [0 0.5 0]);
plot(tnum, num(:, 1), '--', 'linewidth', 4, 'color', [0.7 0 1]);
plot(tnum, num(:, 2), '--',  'linewidth', 4, 'color', [0.7 1 0]);
