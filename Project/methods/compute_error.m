function error = L2_error(K,N,x,D,lambda,R,v)
FinalTime = 0.1;
x1 = -0.2;
x2 = -x1;
C0_peak = 1;
D = 0.01;

C0 = C0fun(x,C0_peak,x1,x2);
[C] = ADR1D_withq(C0,FinalTime,D);
error_vec = C-true_sol(x,FinalTime,v,D,lambda,R,C0_peak,x2);
error = max(error_vec, [], 'all');

function C0 = C0fun(x,magnitude,x1,x2)
    C0 = zeros(size(x));
    ids = (x>x1) & (x<x2);
    C0(ids) = magnitude;
end

function u = true_sol(x,t,v,D,lambda,R,C0_peak,x2)
    u = exp(-lambda*R*t)*(C0_peak/2)*(erf((x-v*t+x2)/sqrt(4*D*t))-erf((x-v*t-x2)/sqrt(4*D*t)));
end

end

