clear
clc
close all

data = dlmread('out.dat');
plot(data(:, 1), data(:, 2));
xlabel('y');
ylabel('u');

plot(data(:, 1), data(:, 3));
xlabel('y');
ylabel('T');
