clear
clc
close all

y0 = 0;
y1 = 1;

k = 1.3;
dy = (y1-y0);
a = 0.1;
y = linspace(y0, y1, 500);

f = tanh(k*((y+a*dy)-y0)/dy) + tanh(k*(y1-(y-a*dy))/dy) - 1;


alpha0 = k/dy;
alpha1 = -k/dy;
beta0  =  alpha0*(a*dy-y0);
beta1  =  alpha0*(y1+a*dy);
    
fi = log(abs(cosh(alpha0*y +beta0)))/alpha0 + log(abs(cosh(alpha1*y +beta1)))/alpha1 - y;
f0 = log(abs(cosh(alpha0*y0+beta0)))/alpha0 + log(abs(cosh(alpha1*y0+beta1)))/alpha1 - y0;
f1 = log(abs(cosh(alpha0*y1+beta0)))/alpha0 + log(abs(cosh(alpha1*y1+beta1)))/alpha1 - y1;
    
plot(y, y0 + dy*(fi-f0)/(f1-f0))

yy = 0;
ff0 = y0 + dy*((log(abs(cosh(alpha0*yy+beta0)))/alpha0 + log(abs(cosh(alpha1*yy+beta1)))/alpha1 - yy)-f0)/(f1-f0);
yy = 0.25;
ff1 = y0 + dy*((log(abs(cosh(alpha0*yy+beta0)))/alpha0 + log(abs(cosh(alpha1*yy+beta1)))/alpha1 - yy)-f0)/(f1-f0);
yy = 0.5;
ff2 = y0 + dy*((log(abs(cosh(alpha0*yy+beta0)))/alpha0 + log(abs(cosh(alpha1*yy+beta1)))/alpha1 - yy)-f0)/(f1-f0);
yy = 0.75;
ff3 = y0 + dy*((log(abs(cosh(alpha0*yy+beta0)))/alpha0 + log(abs(cosh(alpha1*yy+beta1)))/alpha1 - yy)-f0)/(f1-f0);
yy = 1;
ff4 = y0 + dy*((log(abs(cosh(alpha0*yy+beta0)))/alpha0 + log(abs(cosh(alpha1*yy+beta1)))/alpha1 - yy)-f0)/(f1-f0);
[ff0 ff1 ff2 ff3 ff4]