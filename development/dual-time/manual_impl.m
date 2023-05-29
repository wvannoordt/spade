clear
clc
close all
format long

q0 =  5 ;
q1 =  5.07845909572561 ;
q2 =  5.1564344650358 ;
t0 =  0.157079632675 ;
t1 =  0.2356194490125 ;

dt = t1 - t0;
dtau = dt/10.0;
w0 =  ftr(q0);
w1 =  ftr(q1);
w2 =  ftr(q2);

q = q2;
q_ana = 5.0 + sin(t1);

%need to converge an RK iteration
eps = 1000.0;
max_its = 100000;
its = 0;
tau = 0.0;

[rhs_ini, gt] = rk_rhs(q, w0, w1, w2, t1, tau, dt);

while ((abs(eps)>1e-8) && (its < max_its))
    
    [rhs1 gt1] = rk_rhs(q, w0, w1, w2, t1, tau, dt);
    qnew = itr(ftr(q) + 0.5*dtau*rhs1);
    [rhs2 gt2] = rk_rhs(qnew, w0, w1, w2, t1, tau+0.5*dtau, dt);
    dtau
    [rhs1 rhs2]
    [q, qnew]
    error('a');
    eps = rhs2;
    q = ftr(q);
    q = q + dtau*eps;
    q = itr(q);
    its = its + 1;
    tau = tau + dtau;
end



function [qout] = ftr(q)
    qout = q*q;
end

function [qout] = itr(q)
    qout = sqrt(q);
end

function [rhsout] = rhs(q, t)
    rhsout = 2*q*cos(t);
end

function [rk_rhs_out, gt] = rk_rhs(q, w0, w1, w2, t, tau, dt)
    gt = -(18.0/11)*w2 + (9.0/11)*w1 - (2.0/11.0)*w0;
    rk_rhs_out = -ftr(q) - gt + (6.0/11)*dt*rhs(q,t);
end