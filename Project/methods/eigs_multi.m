close all
clearvars
addpath(genpath(pwd))
Globals1D;
upwind = false;
v = 1;
lambda = 1;
R = 1;
N = 1;
K = 20;
D = 0.01;
alp = 1; % 1 is central flux and 0 is upwind flux
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

semilogx(real(eig(LN)),imag(eig(LN)),'o')






