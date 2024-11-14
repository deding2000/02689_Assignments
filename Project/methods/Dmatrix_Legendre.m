function D = Dmatrix_Legendre(N,xjs,a,b)
    dV = GradJacobiP(xjs,0,0,N,true);
    V = JacobiP(xjs,0,0,N,true)';
    Vi = inv(V);
    D = dV*Vi * 2/(b-a); % shifted interval
end