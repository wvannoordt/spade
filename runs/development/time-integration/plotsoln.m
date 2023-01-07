clear
clc
close all

data = dlmread('soln.dat');
figure
subplot(2, 1, 1)
plot(data(:,1), data(:,2), 'linewidth', 4)
hold on
plot(data(:,1), data(:,3), '--', 'linewidth', 3)
legend("numerical", "analytical");
subplot(2, 1, 2)
plot(data(:,1), data(:,2)-data(:,3), 'linewidth', 3)
legend("error");
