clear
clc
close all

data = dlmread('soln.dat');
figure
plot(data(:,1), data(:,2)+(data(1,3)-data(1,2)))
hold on
plot(data(:,1), data(:,3))