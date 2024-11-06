import numpy as np
import scipy as sp
from scipy.special import gamma
import matplotlib.pyplot as plt

def diff_h(N,j,x,xjs):
    dx = x-xjs[j]
    dh = np.cos(N*dx/2)/np.tan(dx/2)/2
    dh += np.sin(N*dx/2)*(-(1+1/np.tan(dx/2)**2)/2)/N
    return dh

def Dmatrix_Fourier(N,xjs):
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
        if isinstance(x, (int, float)):
            lx = 1
        else:
            lx = len(x) 
        P = np.empty((n+1,lx))
        P[0,:] = JacobiP(x, alpha, beta, 0, matrix=False)
        if n>0:
            P[1] = JacobiP(x, alpha, beta, 1, matrix=False)
        if n>1:
            for p in range(2,n+1):
                P[p] = (( a_coef(n=p-1,m=0) + x )*P[p-1] - a_coef(n=p-1,m=-1)*P[p-2])/a_coef(n=p-1,m=1)
        
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
    # if alpha**2 - beta**2 > 10 * np.finfo(float).eps: # added by me - otherwise diagonal is inf
    #J += 
    J += np.diag(-1/2 * (alpha**2 - beta**2) *(1/ (h1 + 2)) *(1/ h1)) + np.diag(np.array(2 / (h1[:N] + 2) * np.sqrt((np.arange(1, N + 1) * 
    (np.arange(1, N + 1) + alpha + beta) * 
    (np.arange(1, N + 1) + alpha) * 
    ((np.arange(1, N + 1) + beta) *  
    1/((h1[:N] + 1)) *(1/ (h1[:N] + 3)))))), 1)

    if alpha + beta < 10 * np.finfo(float).eps:
        J[0, 0] = 0.0
    J += J.T
    #print(J)
    # Compute quadrature by eigenvalue solve
    D, V = np.linalg.eig(J)
    idx = D.argsort()#[::]   
    x = D[idx]
    V = V[:,idx]
    if (alpha == 0) and  (beta == 0):
            w = (V[0, :] ** 2) * 2
    else:
        w = (np.square(V[0, :].T)) * (2 ** (alpha + beta + 1)) / (alpha + beta + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(alpha + beta + 1)
    
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

def GradJacobiP(x,alpha,beta,n,matrix=False):
    if matrix:
        JP = JacobiP(x,alpha+1,beta+1,n-1,matrix=True) 
        vec = (alpha+beta+np.arange(1,n+1)+1)/2
        dV = JP * vec[:,None]
        dV = np.block([[np.zeros(np.shape(dV)[1])],[dV]])
        return  dV.T
    else:
        return (alpha+beta+n+1)/2 * JacobiP(x,alpha+1,beta+1,n-1)
    
def LTM_2ord(X, N, a, b, ffun):
    # Legendre Tau Method
    # Lu = au'' + bu' = 1
    # homogenous boundary conditions

    x = JacobiGL(0,0,N)

    A = np.zeros((N+1,N+1))
    for i in range(2,N):
        A[i,[i-1,i,i+1]] = [ b/a/(2*(i)-1) , 1 , -b/a/(2*(i)+3) ]
    A[N,[N-1,N]] = [b/a/(2*N-1) , 1]
    #nvec = 2/(2*np.arange(N+1) + 1)
    A[0] = (JacobiP(-1,0,0,N,matrix=True)).T#/nvec
    A[1] = (JacobiP(1,0,0,N,matrix=True)).T#/nvec

    # constant f
    f = np.zeros(len(x)+2)
    f[0] = 1
    g = np.empty(N+1)
    g[[0,1]] = [0,0]
    for n in range(2,N+1):
        g[n] = ( f[n-2]/(2*n-3)/(2*n-1) + f[n]*(1/(2*n+3) - 1/(2*n-1)) - f[n+2]/(2*n+5)/(2*n+3) )/a

    u = np.linalg.solve(A,g)

    P = (JacobiP(X,0,0,N,matrix=True)).T#/nvec

    return P@u

def Dmatrix_Legendre(N,xjs,a,b):
        dV = GradJacobiP(xjs,0,0,N,matrix=True)
        V = JacobiP(xjs,0,0,N,matrix=True).T
        Vi = np.linalg.inv(V)
        D = dV@Vi * 2/(b-a) # shifted interval
        return D

def SBM( P, ua, ub, ffun, method, tau = 1, d = 0):

    r = JacobiGL(0,0,P)
    f = ffun(r)
    nvec = np.sqrt(2/(2*np.arange(P+1)+1))
    V = JacobiP(r,0,0,P,matrix=True).T
    Vn = V/nvec
    Vr = GradJacobiP(r,0,0,P,matrix=True)
    M = np.linalg.inv(Vn@Vn.T)
    Vi = np.linalg.inv(V)
    Vin = np.linalg.inv(Vn)
    D = Vr@Vi

    A = D.T @ M @ D 
    b = M @ f

    match method:
        case 'strong_CBM':
            A[0] = 0
            A[0,0] = 1
            b[0] = ua
            A[-1] = 0
            A[-1,-1] = 1
            b[-1] = ub
        case 'weak_CBM':
            A[0] += D[0]
            A[-1] -= D[-1]
            A[0,0] += tau
            A[-1,-1] += tau
            b[0] += tau*ua
            b[-1] += tau*ub
        case 'weak_SBM':
            Phi = JacobiP(np.array([-1-d,1+d]),0,0,P,matrix=True).T
            H = (Phi @ Vi)

            A[0] += D[0] + tau*H[0]
            A[-1] += -D[-1] + tau*H[-1]
            b[0] += tau*ua
            b[-1] += tau*ub
        case _:
            print("No BC method.")
            return
    
    u = np.linalg.solve(A,b)
    condA = np.linalg.cond(A)
    
    return u, r, condA
    



    






