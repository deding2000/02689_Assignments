close all
clearvars
addpath(genpath(pwd))
Globals1D;
v = 1;
lambda = 1;
upwind = false;
R = 1;
Ns = [1,4];
Ks = [10:3:80];
figure
for upw = 1:2
    upwind = ~upwind;
for k=1:length(Ns)
    timesteps = zeros(1,length(Ks));
for j = 1:length(Ks)
N = Ns(k);
K = Ks(j);
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
max_eig = max(abs(eig(LN)));
timesteps(j) = 3/max_eig; 

end

if upwind
    text = 'upwind flux';
else
    text = 'central flux';
end
semilogy(Ks,timesteps,'DisplayName', sprintf('N=%d, %s', N, text))
hold on
end
end
semilogy(Ks,Ks.^(-1.5),'--','DisplayName','1/(K*sqrt(K))')
xlabel('K')
ylabel('\Delta t')
legend('show')


%plot(real(eig(LN)),imag(eig(LN)),'o')









