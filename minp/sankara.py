#
import cvxopt
import numpy as np
from cvxopt import solvers, matrix
#
from topopt_K import K
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
    [Ks,Ds,fs,us,n] = K(nelx,nely,rmin,penal)
#
    u = np.zeros(n,dtype=float)
#
    sol=minimize(func,u,args=(Ks,Ds,fs), \
        jac=grad,method='L-BFGS-B',options={'disp':True, 'gtol': 1e-16, 'ftol': 1e-16, \
        'maxfun':100000, 'maxiter':100000, 'maxcor': n})
#
    us=np.array(us).flatten()
    print(sol.x-us)
#
