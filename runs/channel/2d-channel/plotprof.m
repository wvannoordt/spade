clear
clc
close all

delta = 1;
re_tau = 180;
u_tau = 0.0318736;
data = load('prof.dat');
figure
plot(data(:,1), data(:,2), 'linewidth', 3)
hold on
y = linspace(-1,1,300);
plot(y, 0.5*re_tau*u_tau*(1-y.^2/delta^2), '--', 'linewidth', 3);