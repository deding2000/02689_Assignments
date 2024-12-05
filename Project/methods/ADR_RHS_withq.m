function [rhsu] = ADR_RHS_withq(u,D,upwind)
% function [rhsu] = AdvecRHS1D_withq(u,time)
% Purpose : Evaluate RHS flux in 1D ADR equation with first order rewrite

Globals1D;
% form field differences at faces
du = zeros(Nfp*Nfaces,K); du(:) = u(vmapM)-u(vmapP);

% impose boundary condition at x=-1,1
uin = 0;
du(mapI) = (u(vmapI)-uin); %*2.0
uout=0;
du(mapO) = (u(vmapO)-uout); %*2.0

% Compute q with central flux
q = sqrt(D)*(rx.*(Dr*u) - LIFT*(Fscale.*(nx.*du/2.0)));
dq = zeros(Nfaces,K); dq(:) = (q(vmapM)-q(vmapP))/2.0;

% impose boundary condition - Dirichlet conditions
dq(mapI) = 0.0; dq(mapO) = 0.0;

% Flux for advection term
du2 = zeros(Nfp*Nfaces,K); 
if not(upwind)
% Central flux
du2(:) = (u(vmapM)-u(vmapP))/2.0;
du2(mapI)=(u(vmapI)-uin); du2(mapO)=(u(vmapO)-uout);
flux = nx.*(du2*v - sqrt(D)*dq) - 1/2.0.*du;
end

% Upwind flux
if upwind
du2(:) = (u(vmapM)-u(vmapP)).*(v*nx(:)-(1-0)*abs(v*nx(:)))/2;
du2 (mapI) = (u(vmapI)- uin ).*(v*nx(mapI)-(1-0)*abs(v*nx(mapI)))/2;
flux = -nx.*(sqrt(D)*dq) - du/2.0 + du2;
end

% Local derivatives
dfdr = Dr*(u*v - sqrt(D)*q);

% compute right hand sides of the semi-discrete PDE
rhsu = -(rx.*dfdr - LIFT*(Fscale.*flux))-lambda*R*u;
return