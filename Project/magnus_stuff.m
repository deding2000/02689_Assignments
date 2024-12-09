%%
close all
clearvars
addpath(genpath(pwd))

%%
N = 40;
T = 0.001;
D = 1;
v = 5;
R  = 1;
lam = 1;
gfun = @(t) 0;
ffun = @(t) 0;




[t,C,x] = ADRsolver1D(N,T,D,v,lam,R,gfun,ffun,false);



%%
figure
plot(x,C(end,:))

%% 
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


Ns = [2]*5;
maxe = zeros(size(Ns));
for i = 1:length(Ns)
    N = Ns(i);
    [~,~,~,sys] = ADRsolver1D(N,3,D,v,lambda,R,@(x)0,@(x)0,C0);
    maxe(i) = max(abs(real(eig(sys))));
end

figure
semilogy(Ns,maxe)
hold on 
plot(Ns,Ns.^(2.5),'--')
legend('Max. eigenvalue','reference, N^{5/2}')
xlabel('N')
ylabel('\lambda')

%% 
close all
clearvars
addpath(genpath(pwd))
Globals1D;
% global N Nfp Np K
% global r x VX
% global Dr LIFT
% global nx Fx Fscale
% global vmapM vmapP vmapB mapB Fmask
% global vmapI vmapO mapI mapO
% global rx J
% global rk4a rk4b rk4c
% global Nfaces EToE EToF
% global V invV
% global NODETOL
% % Parameters of model
% global D v lambda R

D = 0;
v = 1;
lambda = 1;
R = 1;

N = 3;
K = 10;
StartUp1D;

FinalTime = 0.01;
C0 = C0fun(x,1,-0.2,0.2);

figure(1)
plot(x,C0)
hold on

[C] = ADR1D(C0,FinalTime);

figure(1)
plot(x,C)




function C0 = C0fun(x,magnitude,x1,x2)
    C0 = zeros(size(x));
    ids = (x>x1) & (x<x2);
    C0(ids) = magnitude;
end

function u = true_sol(x,t,v,D,lambda,R,C0_peak,x2)
    u = exp(-lambda*R*t)*(C0_peak/2)*(erf((x-v*t+x2)/sqrt(4*D*t))-erf((x-v*t-x2)/sqrt(4*D*t)));
end


