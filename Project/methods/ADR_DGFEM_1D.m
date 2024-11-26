function [X] = ADR_DGFEM_1D(K,N,T,eq_param,fluxA,fluxD)

VX = linspace(-1,1,K+1);
h = 2/K; %diff(VX);

% unpack parameters
D = eq_param.D;
v = eq_param.v;
lam = eq_param.lam;
R = eq_param.R;
gfun = eq_param.gfun;
ffun = eq_param.ffun;
C0 = eq_param.C0;

xr = JacobiGL(0,0,N);
V = JacobiP(xr,0,0,N+1,true)';
nvec = sqrt( 2./(2.*(0:(N+1))+1));
V = V ./ nvec;
M = inv(V*V');
Mk = h/2*M;
Mki = inv(Mk);
Vx = GradJacobiP(xr,0,0,N+1,true)';
Vx = Vx./nvec;
Vi = inv(V);
Dx = Vx*Vi;
S = M*Dx;
S2 = S*Dx;


k_int = D*S2 - v*S - lam*R*M;
RHS = zeros(K*(N+1));






