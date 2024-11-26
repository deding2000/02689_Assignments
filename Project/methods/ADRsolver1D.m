function [t,C,x] = ADRsolver1D(N,T,D,v,lam,R,gfun,ffun,C0)
% dC/dt = D * d^2C/dx^2 - v * dC/dx - lam*R * C
%   C(x,0) = C0(x)
%   C(-1,t) = gfun(t)
%   (dC/dx)(-1,t) = ffun(t)

[x,w] = JacobiGLw(0,0,N);
%C0 = %arrayfun(C0fun,x);
Dm = Dmatrix_Legendre(N,x,-1,1);
W = diag(w);
Wi = diag(1/w);
H = zeros(size(W));
H(end,end) = 1;
e0 = zeros(N+1,1);
e0(1) = 1;

dCdt = @( t, C ) ( Wi*(D*Dm - v*H) - D*Dm'*Dm + v*Dm' - lam*R*eye(N+1) )*C - D*Wi*e0*ffun(t) + v*Wi*e0*gfun(t);
tspan = [0 T];

if ~C0
    C0 = zeros(size(x));
    C0(abs(x)<0.3) = 1;
end

[t,C] = ode45(dCdt,tspan,C0);


return




