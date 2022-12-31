#
import cvxopt
import numpy as np
from cvxopt import solvers, matrix
from scipy.sparse import hstack
from scipy.sparse.linalg import eigs
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
    dofs=np.arange(2*(nelx+1)*(nely+1))
    ndof = 2*(nelx+1)*(nely+1)
    fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
    x=np.ones(nelx*nely,dtype=float)
    n=ndof-len(fixed)
    u = np.zeros(n,dtype=float)
#
    [Hs,Ds,Gs] = K(nelx,nely,rmin,penal,x,u,fixed)
#
    ones = np.append(np.ones(n+nelx*nely),-np.ones(n+nelx*nely))
    J = cvxopt.spmatrix(ones, range(2*n+2*nelx*nely), np.append(range(n+nelx*nely),range(n+nelx*nely)))
    ones = np.append(np.ones(n+nelx*nely),np.ones(n+nelx*nely))
    h = cvxopt.matrix(ones)
#
    sol=solvers.qp(Hs, Gs, J, h)
#
    u[:]=sol['x'][:n]
    x[:]=sol['x'][n:]
#
    [Hs,Ds,Gs] = K(nelx,nely,rmin,penal,x,u,fixed)
#
