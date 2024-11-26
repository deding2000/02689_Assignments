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


