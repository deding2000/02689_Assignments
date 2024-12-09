function D = Dmatrix_Legendre(N,xjs,a,b)
    dV = GradVandermonde1D(N,xjs);
    %V = JacobiP(xjs,0,0,N)';
    V = Vandermonde1D(N,xjs);
    Vi = inv(V);
    D = dV*Vi * 2/(b-a); % shifted interval
end