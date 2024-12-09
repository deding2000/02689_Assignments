function [rhsu] = ADR_RHS_withq(u,D,upwind)
% function [rhsu] = AdvecRHS1D_withq(u,time)
% Purpose : Evaluate RHS in 1D ADR equation with first order rewrite
Globals1D; 
% form field differences at faces
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);
% impose boundary conditions at x=-1,1
uin = 0;
du(mapI) = (u(vmapI)-uin); 
uout=0;
du(mapO) = (u(vmapO)-uout);

% Compute q with central flux
q = sqrt(D)*(rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du/2.0)));
dq = zeros(Nfaces,K); dq(:) = (q(vmapM)-q(vmapP))/2.0;

% impose boundary condition - Dirichlet conditions
dq(mapI) = 0.0; dq(mapO) = 0.0;

% Central or upwind flux for advection term
if upwind
alpha = 0;
else
alpha = 1;
end

du2 = zeros(Nfp*Nfaces,K); 
du2(:) = (u(vmapM)-u(vmapP)).*(v*nx(:)-(1-alpha)*abs(v*nx(:)))/2;
du2(mapI) = (u(vmapI)- uin ).*(v*nx(mapI)-(1-alpha)*abs(v*nx(mapI)))/2;
du2 (mapO) = 0;

% Total flux
flux = -nx.*(sqrt(D)*dq) - du/2.0 + du2;

% Local derivatives
dfdr = Dr*(v*u - sqrt(D)*q);

% compute right hand sides of the semi-discrete PDE
rhsu = -(rx.*dfdr - LIFT*(Fscale.*flux))-lambda*R*u;
return