%close all
clearvars
addpath(genpath(pwd))
Globals1D;
upwind = false; % choose upwind flux or not for the advection term
v = 1;
lambda = 1;
R = 1;
N = 3;
K = 10;
D = 0.01;
StartUp1D;
FinalTime = 0.1;
BIGN = K*(N+1);
LN = zeros(BIGN,BIGN);
for i = 1:BIGN
ev = zeros(N+1,K);
ev(i) = 1;
[rhsu] = ADR_RHS_withq(ev,D,upwind);
LN(:,i) = rhsu(:);
end
figure
plot(real(eig(LN)),imag(eig(LN)),'o',color="r")
hold on

clearvars
addpath(genpath(pwd))
Globals1D;
upwind = true; % choose upwind flux or not for the advection term
v = 1;
lambda = 1;
R = 1;
N = 3;
K = 10;
D = 0.01;
StartUp1D;
FinalTime = 0.1;
BIGN = K*(N+1);
LN = zeros(BIGN,BIGN);
for i = 1:BIGN
ev = zeros(N+1,K);
ev(i) = 1;
[rhsu] = ADR_RHS_withq(ev,D,upwind);
LN(:,i) = rhsu(:);
end
plot(real(eig(LN)),imag(eig(LN)),'o',color="b")
xlabel('$Re(\lambda)$','Interpreter','latex') 
ylabel('$Im(\lambda)$','Interpreter','latex')
legend({'central','upwind'},'Location','northeast')






