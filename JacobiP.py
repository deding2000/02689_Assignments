import numpy as np

def JacobiP(x,alpha,beta,n):
    if n == 0:
        return np.ones(len(x))
    elif n == 1:
        return 1/2*(alpha-beta+(alpha+beta+2)*x)
    else:
        a1 = (2*(n+alpha)*(n+beta))/((2*n+alpha+beta+1)*(2*n+alpha+beta))
        ann = (alpha**2-beta**2)/((2*n+alpha+beta+2)*(2*n+alpha+beta))
        a2 = (2*(n+1)*(n+alpha+beta+1))/((2*n+alpha+beta+2)*(2*n+alpha+beta+1))
        return ((ann+x)*JacobiP(x,alpha,beta,n-1)-a1*JacobiP(x,alpha,beta,n-2))/a2