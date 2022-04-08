clear
clc
close all

xrq = [0.0520833 1.55208 -0.447917]

xyz = [xrq(1) xrq(2)*cos(xrq(3)) xrq(2)*sin(xrq(3))]

m = [1 0 0; 0 cos(xrq(3)) -xrq(2)*sin(xrq(3)); 0 sin(xrq(3)) xrq(2)*cos(xrq(3))]

J = det(m^-1)