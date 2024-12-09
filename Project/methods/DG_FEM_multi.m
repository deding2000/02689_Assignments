%% 
% Solve ADR equation in 1D using the DG-FEM method with multiple domains
close all
clearvars
addpath(genpath(pwd))
Globals1D;
upwind = true; %choose upwind flux or not for the advection term
v = 1;
lambda = 1;
R = 1;
N = 1;
K = 30;
D = 0.01;
StartUp1D;
FinalTime = 0.1;
x1 = -0.2;
x2 = -x1;
C0_peak = 1;
C0 = C0fun(x,C0_peak,x1,x2);
[C] = ADR1D_withq(C0,FinalTime,D,upwind);
% Print some errors
disp(norm(C-true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2)))
disp(max(C-true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2), [], 'all'))


% plot true sol and approxmations
figure(1)
xx = linspace(-1,1);
true_C = true_sol(xx,0,v,D,lambda,R,C0_peak,x2);
plot(xx,true_C,'Color','black')
hold on
true_C = true_sol(xx,FinalTime,v,D,lambda,R,C0_peak,x2);
plot(xx,true_C,'Color','black','LineStyle','--')
hold on
plot(x,C,'Color','blue')

hold on
FinalTime = 0.2;
true_C = true_sol(xx,FinalTime,v,D,lambda,R,C0_peak,x2);
plot(xx,true_C,'black','LineStyle','--')
hold on
[C2] = ADR1D_withq(C0,FinalTime,D,upwind);
plot(x,C2,'blue')
hold on
FinalTime = 0.3;
true_C = true_sol(xx,FinalTime,v,D,lambda,R,C0_peak,x2);
plot(xx,true_C,'black','LineStyle','--')
hold on
[C3] = ADR1D_withq(C0,FinalTime,D,upwind);
plot(x,C3,'blue')
hold off
xlabel('$x$','Interpreter','latex') 
ylabel('$C$','Interpreter','latex')
legend({'$C(x,0)$','$C(x,t)$','$C_N(x,t)$'},'Location','northeast','Interpreter','latex')

function C0 = C0fun(x,magnitude,x1,x2)
    C0 = zeros(size(x));
    ids = (x>x1) & (x<x2);
    C0(ids) = magnitude;
end

function u = true_sol(x,t,v,D,lambda,R,C0_peak,x2)
    u = exp(-lambda*R*t)*(C0_peak/2)*(erf((x-v*t+x2)/sqrt(4*D*t))-erf((x-v*t-x2)/sqrt(4*D*t)));
end







