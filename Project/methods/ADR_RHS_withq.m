function [rhsu] = ADR_RHS_withq(u,D)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose : Evaluate RHS flux in 1D ADR equation with first order rewrite
% COMPARE WITH page 256 in yellow book
Globals1D;
% form field differences at faces
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);

% impose boundary condition at x=-1,1
uin = 0;
du(mapI) = 2.0*(u(vmapI)-uin);
uout=0;
du(mapO) = 2.0*(u(vmapO)-uout);

% Compute q and jumps
q = sqrt(D)*(rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du/4.0)));
dq = zeros(Nfaces,K); dq(:) = (q(vmapM)-q(vmapP))/2.0;

% impose boundary condition - Dirichlet conditions
dq(mapI) = 0.0; dq(mapO) = 0.0;

% Evaluate nonlinear flux (maybe Lax-Friedrich flux - we should use something
% else?)
du2 = zeros(Nfp*Nfaces,K); du2(:) = (u(vmapM)-u(vmapP))/2.0;

% impose boundary condition
du2(mapI)=(u(vmapI).^1-uin.^1); du2(mapO)=(u(vmapO).^1-uout.^1);

% Compute flux
maxvel = max(max(abs(u)));

% penalty scaling -- See Chapter 7.2 ????
%tau = .25*reshape(N*N./max(2*J(vmapP),2*J(vmapM)), Nfp*Nfaces, K);
tau=0;

% flux term
flux = nx.*(du2*v - sqrt(D)*dq) - maxvel/2.0.*du ...
- sqrt(D)*tau.*du;

% local derivatives of field
dfdr = Dr*(u*v - sqrt(D)*q);

% compute right hand sides of the semi-discrete PDE
rhsu = -(rx.*dfdr - LIFT*(Fscale.*flux))-lambda*R*u;
return