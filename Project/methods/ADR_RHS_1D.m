function [rhsu] = ADR_RHS_1D(u,time, v)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose : Evaluate RHS flux in 1D advection
Globals1D;
% form field differences at faces
alpha=1;
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(v*nx(:)-(1-alpha)*abs(v*nx(:)))/2;
% impose boundary condition at x=0
uin = 0*time; %-sin(a*time);
du (mapI) = (u(vmapI)- uin ).*(v*nx(mapI)-(1-alpha)*abs(v*nx(mapI)))/2;
du (mapO) = 0;
% compute right hand sides of the semi-discrete PDE
advec = -v*rx.*(Dr*u) + LIFT*(Fscale.*(du));
diff = D*rx.*(Dr*Dr*u) + LIFT*(Fscale.*(du));
rhsu = advec;
return