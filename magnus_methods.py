import numpy as np
import scipy as sp
from scipy.special import gamma
import matplotlib.pyplot as plt
from methods import *

def massGL(N):
    x = JacobiGL(0,0,N)
    nvec = np.sqrt(2/(2*np.arange(N+1)+1))
    V = JacobiP(x,0,0,N,matrix=True).T/nvec
    return np.linalg.inv(V@V.T)

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
