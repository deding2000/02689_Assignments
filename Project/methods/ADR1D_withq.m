function [u] = ADR1D_withq(u, FinalTime,D,upwind)
% function [u] = Advec1D(u, FinalTime)
% Purpose : Integrate 1D ADR equation until FinalTime starting with
% initial the condition, u
Globals1D;
time = 0;
% Runge-Kutta residual storage
resu = zeros(Np,K);
% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.25; % 0.75 before
dt = CFL*(xmin)^2; 
dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps;
% outer time step loop
for tstep=1:Nsteps
    for INTRK = 1:5
        timelocal = time + rk4c(INTRK)*dt;
        %[rhsu] = ADR_RHS_1D(u, timelocal);
        [rhsu] = ADR_RHS_withq(u,D,upwind);
        resu = rk4a(INTRK)*resu + dt*rhsu;
        u = u+rk4b(INTRK)*resu;
    end
    % Increment time
    time = time+dt;
end
return