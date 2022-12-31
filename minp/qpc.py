#
import cvxopt
import numpy as np
from cvxopt import solvers, matrix
#
from topopt_H import K
#
from scipy.optimize import minimize
#
def func(u,Ks,Ds,fs):
#
    us = cvxopt.matrix(u)
    R = Ks*us - fs
#
    return R.T * Ds * R
#
def grad(u,Ks,Ds,fs):
#
    us = cvxopt.matrix(u)
    R = Ks*us - fs
#
    return 2.* R.T * Ds * Ks 
#
def hess(u,Ks,Ds,fs):
#
    us = cvxopt.matrix(u)
    R = Ks*us - fs
#
    return 2.* Ks.T * Ds * Ks
#
if __name__ == "__main__":
#
    nelx=18
    nely=6
    rmin=2.
    penal=3.0
#
    [Ks,Ds,fs,u0,n] = K(nelx,nely,rmin,penal)
#
    u = np.zeros(n,dtype=float)
#
    ones = np.append(np.ones(n),-np.ones(n))
    J = cvxopt.spmatrix(ones, range(2*n), np.append(range(n),range(n)))
    ones = np.append(np.ones(n),np.ones(n))
    h = cvxopt.matrix(ones*10.)
#
    for k in range(100):   
#
        us = cvxopt.matrix(u)
        [Ks,Ds,fs,_,n] = K(nelx,nely,rmin,penal)
        sol=solvers.qp(Ds, Ks*us-fs, J, h)
        u[:] = u + np.array(sol['x']).flatten()
#   
        if np.linalg.norm(sol['x']) < 1e-6:
            break
#
    print(u[0])
#
