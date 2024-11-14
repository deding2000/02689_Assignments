function JP = GradJacobiP(x,alpha,beta,n,matrix)
    if matrix
        JP = JacobiP(x,alpha+1,beta+1,n-1,true);
        vec = (alpha+beta+(1:n)+1)/2;
        dV = JP .* vec';
        dV = [zeros(1, size(dV,2)); dV];
        JP = dV';
    else
        JP = (alpha+beta+n+1)/2 * JacobiP(x,alpha+1,beta+1,n-1,false);
    end
end