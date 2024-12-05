%% 
close all
clearvars
addpath(genpath(pwd))
Globals1D;
v = 1;
lambda = 1;
R = 1;
N = 2;
K = 40;
D = 0.01;
alp = 1; % 1 is central flux and 0 is upwind flux
StartUp1D;
FinalTime = 0.1;
x1 = -0.2;
x2 = -x1;
C0_peak = 1;
C0 = C0fun(x,C0_peak,x1,x2);
[C] = ADR1D_withq(C0,FinalTime,D);

figure(1)
xx = linspace(-1,1);
for t = linspace(0,FinalTime,2)
    true_C = true_sol(xx,t,v,D,lambda,R,C0_peak,x2);
    plot(xx,true_C)
    hold on
end

%figure(2)
%plot(xx,true_C)
%plot(x,C0)
hold on
plot(x,C,'black')
xlabel('$x$','Interpreter','latex') 
ylabel('$C$','Interpreter','latex')
legend({'$C(x,0)$','$C(x,0.1)$','$C_N(x,0.1)$'},'Location','northeast','Interpreter','latex')

function C0 = C0fun(x,magnitude,x1,x2)
    C0 = zeros(size(x));
    ids = (x>x1) & (x<x2);
    C0(ids) = magnitude;
end

function u = true_sol(x,t,v,D,lambda,R,C0_peak,x2)
    u = exp(-lambda*R*t)*(C0_peak/2)*(erf((x-v*t+x2)/sqrt(4*D*t))-erf((x-v*t-x2)/sqrt(4*D*t)));
end

% Convergence test
%KS = [10,20,30];
% errors = [];
% true_C = true_sol(x,t,v,D,lambda,R,C0_peak,x2);
% for KK = KS
%     K = KK;
%     [C] = ADR1D(C0,FinalTime,D);
%     errors = [errors,norm(C-true_C)];
% end

% figure(3)
% plot(KS,errors)
disp(norm(C-true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2)))
disp(max(C-true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2), [], 'all'))





