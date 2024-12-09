function [rhsu] = ADR_RHS_1D(u,time)
% function [rhsu] = AdvecRHS1D(u,time)
% Purpose : Evaluate RHS flux in 1D advection
Globals1D;
% form field differences at faces
du = zeros(Nfp*Nfaces,K);
du(:) = (u(vmapM)-u(vmapP)).*(v*nx(:)-(1-alp)*abs(v*nx(:)))/2;
du2 = zeros(Nfp*Nfaces,K);
du2(:) = (u(vmapM)-u(vmapP)).*(D*nx(:));%-(1-alp)*abs(v*nx(:)))/2;
% impose boundary condition at x=0
uin = 0*time; %-sin(a*time);
du (mapI) = (u(vmapI)- uin ).*(v*nx(mapI)-(1-alp)*abs(v*nx(mapI)))/2;
du (mapO) = 0;
% compute right hand sides of the semi-discrete PDE
advec = -v*rx.*(Dr*u) + LIFT*(Fscale.*(du));
diff = D*rx.*(Dr*Dr*u) - LIFT*(Fscale.*(du2)); % du2 ??
reac = -lambda*R*u;
rhsu = advec+reac;
return