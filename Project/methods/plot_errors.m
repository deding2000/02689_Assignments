close all
clearvars
addpath(genpath(pwd))
Globals1D;
v = 1;
lambda = 1;
R = 1;

D = 0.01;
alp = 1;

Ns = [1,2,3,4];
KS = [20,30,50,100];
figure
for NN = Ns
N = NN;
error = [];
for KK = KS
K = KK;
% STARTUP
% Purpose : Setup script, building operators, grid, metric and
% connectivity for 1D solver.
% Definition of constants
Globals1D; NODETOL = 1e-10;
Np = N+1; Nfp = 1; Nfaces=2;
% Compute basic Legendre Gauss Lobatto grid
VX = linspace(-1,1,K+1);
EToV = [1:K;2:K+1]';
r = JacobiGL(0,0,N);
% Build reference element matrices
V = Vandermonde1D(N, r); invV = inv(V);
Dr = Dmatrix1D(N, r, V);
% Create surface integral terms
LIFT = LIFT1D();
% build coordinates of all the nodes
va = EToV(:,1)'; vb = EToV(:,2)';
x = ones(N+1,1)*VX(va) + 0.5*(r+1)*(VX(vb)-VX(va));
% calculate geometric factors
[rx,J] = GeometricFactors1D(x,Dr);
% Compute masks for edge nodes
fmask1 = find( abs(r+1) < NODETOL)';
fmask2 = find( abs(r-1) < NODETOL)';
Fmask = [fmask1;fmask2]';
Fx = x(Fmask(:), :);
% Build surface normals and inverse metric at surface
[nx] = Normals1D();
Fscale = 1./(J(Fmask,:));
% Build connectivity matrix
[EToE, EToF] = Connect1D(EToV);
% Build connectivity maps
[vmapM, vmapP, vmapB, mapB] = BuildMaps1D;
error = [error,compute_error(N,K,x,D,lambda,R,v)];
end

loglog(KS,error)
hold on
end

xlabel('$K$','Interpreter','latex') 
ylabel('$\|C-C_N\|_{\infty}$','Interpreter','latex')
legend({'N=1','N=2','N=3','N=4'},'Location','northeast')
