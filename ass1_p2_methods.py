import numpy as np
import scipy as sp
from scipy.special import gamma

def diff_h(N,j,x,xjs):
    dx = x-xjs[j]
    dh = np.cos(N*dx/2)/np.tan(dx/2)/2
    dh += np.sin(N*dx/2)*(-(1+1/np.tan(dx/2)**2)/2)/N
    return dh

'''def diff_h(N,j,x,xjs):
    dx = x-xjs[j]
    dh = np.cos(N*dx/2)/np.tan(dx/2)
    dh += np.sin(N*dx/2)*(1/2 + np.tan(dx/2)**2)/np.tan(dx/2)**2
    dh /= N
    return dh'''

def Dmatrix(N,xjs):
    D = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
                # if i != 0 and i != (N-1):
                    if i == j:
                        D[i][j] = 0
                    else: 
                        D[i][j] = diff_h(N,j,xjs[i],xjs)
    for i in range(N):
         s = np.sum(D[i,:])
         D[i][i] = -s
    return D

def JacobiP(x, alpha, beta, n, matrix = False):

    def a_coef(n, m):
        # a^(alpha,beta)_(n+m,n)
        if m==-1:
            a = 2*(n+alpha)*(n+beta)/(2*n+alpha+beta+1)/(2*n+alpha+beta)
        elif m==0:
            a = (alpha**2-beta**2)/(2*n+alpha+beta+2)/(2*n+alpha+beta)
        else:#m==1
            a = 2*(n+1)*(n+alpha+beta+1)/(2*n+alpha+beta+2)/(2*n+alpha+beta+1)
        return a

    if not(matrix):
        if n == 0:
            P = np.ones_like(x)
        elif n == 1:
            P = (alpha - beta + (alpha + beta + 2)*x)/2
        else:
            P = ( a_coef(n=n-1,m=0) + x )*JacobiP(x,alpha,beta,n-1) - a_coef(n=n-1,m=-1)*JacobiP(x,alpha,beta,n-2)
            P /= a_coef(n=n-1,m=1)
    else:
        P = np.zeros((n+1,len(x)))
        P[0,:] += JacobiP(x, alpha, beta, 0, matrix=False)
        if n>0:
            P[1] += JacobiP(x, alpha, beta, 1, matrix=False)
        if n>1:
            for m in range(2,n+1):
                P[m] += (( a_coef(n=n-1,m=0) + x )*P[m-1] - a_coef(n=n-1,m=-1)*P[m-2])/a_coef(n=n-1,m=1)
        
    return P

def JacobiGQ(alpha, beta, N):
    # converted code from Allan
    if N == 0:
        x = np.array([-(alpha - beta) / (alpha + beta + 2)])
        w = np.array([2])
        return x, w
    
    # Form symmetric matrix from recurrence.
    J = np.zeros((N + 1, N + 1))
    h1 = 2 * np.arange(N + 1) + alpha + beta
    if alpha**2 - beta**2 > 10 * np.finfo(float).eps: # added by me - otherwise diagonal is inf
        J += np.diag(-1 / (2 * (alpha**2 - beta**2) / (h1 + 2) / h1))
    J += np.diag(np.array(2 / (h1[:-1] + 2) * np.sqrt((np.arange(1, N + 1) * 
        (np.arange(1, N + 1) + alpha + beta) * 
        (np.arange(1, N + 1) + alpha) * 
        (np.arange(1, N + 1) + beta) / 
        (h1[1:] + 1) / (h1[1:] + 3)))), 1)

    if alpha + beta < 10 * np.finfo(float).eps:
        J[0, 0] = 0.0

    J += J.T

    # Compute quadrature by eigenvalue solve
    D, V = np.linalg.eig(J)
    x = np.diag(D)
    w = (V[0, :] ** 2) * (2 ** (alpha + beta + 1)) / (alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 1)
    
    return x, w

def JacobiGL(alpha, beta, N):
    # Purpose: Compute the N'th order Gauss Lobatto quadrature 
    #          points, x, associated with the Jacobi polynomial,
    #          of type (alpha,beta) > -1 ( <> -0.5). 
    # Converted code from Allan

    x = np.zeros(N + 1)
    if N == 1:
        x[0] = -1.0
        x[1] = 1.0
        return x

    xint, _ = JacobiGQ(alpha + 1, beta + 1, N - 2)
    x = np.concatenate(([-1], xint, [1]))
    return x


