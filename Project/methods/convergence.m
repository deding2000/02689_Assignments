%% 
close all
clearvars
addpath(genpath(pwd))
% Convergence test
KS = [10,20,30];
errors = [];
for KK = KS
Globals1D;
K = KK;
disp(K)
v = 1;
lambda = 1;
R = 1;
N = 3;
K = 30;
D = 0.001;
alp = 1; % 1 is central flux and 0 is upwind flux
StartUp1D;
FinalTime = 0.1;
x1 = -0.2;
x2 = -x1;
C0_peak = 1;
C0 = C0fun(x,C0_peak,x1,x2);

[C] = ADR1D(C0,FinalTime,D);
true_C = true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2);
errors = [errors,norm(C-true_C)];
end

figure(1)
plot(KS,errors)

function C0 = C0fun(x,magnitude,x1,x2)
    C0 = zeros(size(x));
    ids = (x>x1) & (x<x2);
    C0(ids) = magnitude;
end

function u = true_sol(x,t,v,D,lambda,R,C0_peak,x2)
    u = exp(-lambda*R*t)*(C0_peak/2)*(erf((x-v*t+x2)/sqrt(4*D*t))-erf((x-v*t-x2)/sqrt(4*D*t)));
end