function P = JacobiP(x, alpha, beta, n, matrix)
    function a = a_coef(n, m)
        if m == -1
            a = 2*(n+alpha)*(n+beta)/(2*n+alpha+beta+1)/(2*n+alpha+beta);
        elseif m == 0
            a = (alpha^2 - beta^2)/(2*n+alpha+beta+2)/(2*n+alpha+beta);
        else
            a = 2*(n+1)*(n+alpha+beta+1)/(2*n+alpha+beta+2)/(2*n+alpha+beta+1);
        end
    end
    
    if ~matrix
        if n == 0
            P = ones(size(x));
        elseif n == 1
            P = (alpha - beta + (alpha + beta + 2).*x)/2;
        else
            P = (a_coef(n-1, 0) + x).*JacobiP(x, alpha, beta, n-1) - a_coef(n-1, -1)*JacobiP(x, alpha, beta, n-2);
            P = P./a_coef(n-1, 1);
        end
    else
        P = zeros(n+1, length(x));
        P(1, :) = JacobiP(x, alpha, beta, 0, false);
        if n > 0
            P(2, :) = JacobiP(x, alpha, beta, 1, false);
        end
        if n > 1
            for p = 2:n
                P(p+1, :) = ((a_coef(p-1, 0) + x).*P(p, :) - a_coef(p-1, -1).*P(p-1, :))./a_coef(p-1, 1);
            end
        end
    end
end