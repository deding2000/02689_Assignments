function [x, w] = JacobiGLw(alpha, beta, N)
% function [x, w] = JacobiGL(alpha, beta, N)
% Purpose: Compute the N'th order Gauss-Lobatto quadrature 
%          points, x, and weights, w, associated with the 
%          Jacobi polynomial of type (alpha, beta) > -1 (<> -0.5).

% Allocate memory for points and weights
x = zeros(N+1,1);
w = zeros(N+1,1);

% Special case: N = 1
if (N == 1)
    x(1) = -1.0; 
    x(2) = 1.0; 
    w(1) = 1.0;
    w(2) = 1.0;
    return;
end

% Compute interior Gauss points and weights using JacobiGQ
[xint, wint] = JacobiGQ(alpha+1, beta+1, N-2);

% Add endpoints -1 and 1
x = [-1; xint; 1];

% Compute weights for the Lobatto points
w(1) = 2 / (N * (N + alpha + beta + 1));
w(N+1) = w(1);

for i = 2:N
    % Weight formula using Jacobi polynomial values
    P0 = JacobiP(x(i), alpha, beta, 0, false); % Evaluate P_0
    w(i) = wint(i-1) / P0^2;            % Normalize using P_0
end

end