%%
close all
clearvars
addpath(genpath(pwd))

%%
N = 60;
T = 0.0001;
D = 10;
v = 5;
R  = 1;
lam = 1;
gfun = @(t) 0;
ffun = @(t) 0;
C0 = C0test(N);


[t,C,x] = ADRsolver1D(N,T,D,v,lam,R,gfun,ffun,C0);

%%
figure
plot(x,C(end,:))

%%
function C = C0test(N)
    x = JacobiGL(0,0,N);
    y = @(x) 10*exp(-0.1/(0.1-x^2));
    C = zeros(size(x));
    C(abs(x)<0.1) = arrayfun(y,x(abs(x)<0.1));
return; end