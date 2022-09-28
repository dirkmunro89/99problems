#
import numpy as np
from scipy.optimize import minimize
#
def con(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    bds=[[1e-6,1e8] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(con_dual,x_d,args=(n,m,x_k,g,dg,d_l,d_u), \
        jac=dcon_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u)
#
#   dq = -np.sum(np.where(dg[0]>0.,dg[0]*(x-x_k),-dg[0]*(1e0/x-1e0/x_k)*(x_k)**2e0))
#   dq = dq + np.sum(np.where(dg[0]>0.,dg[0]*(x_k-x_k),-dg[0]*(1e0/x_k-1e0/x_k)*(x_k)**2e0))
    q_k=g.copy()
    for i in range(m+1):
        q_k[i] = q_k[i] + np.sum(np.where(dg[i]>0.,dg[i]*(x-x_k),-dg[i]*(1e0/x-1e0/x_k)*(x_k)**2e0))
#
#   dq=g[0]-q_k[0]
#   print(g-q_k,dq)
#
    return x,d,q_k
#
# CONLIN: x in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u):
#
    x = np.zeros(n,dtype=np.float64)
#
    tmpp = np.zeros(n,dtype=np.float64)
    tmpn = np.zeros(n,dtype=np.float64)
#
    tmpp[:] = np.where( dg[0] > 0e0, dg[0], np.zeros(n,dtype=float))
    tmpn[:] = np.where( dg[0] <= 0e0, -dg[0], np.zeros(n,dtype=float))
#
    for i in range(m):
        tmpp[:] = tmpp + np.where( dg[i+1] > 0e0, dg[i+1]*x_d[i], np.zeros(n,dtype=float))
        tmpn[:] = tmpn - np.where( dg[i+1] <= 0e0, dg[i+1]*x_d[i], np.zeros(n,dtype=float))
#
    tmpp=np.maximum(tmpp,1e-6)
    tmpn=np.maximum(tmpn,0e0)
#
    x = np.minimum(np.maximum(np.sqrt(tmpn/tmpp)*x_k,dx_l),dx_u)
#
    return x
#
# CONLIN: Dual function value
#
def con_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u)
#
    W = g[0] + np.dot(x_d,g[1:])
    W = W + np.sum(np.where(dg[0]>0.,dg[0]*(x-x_k),-dg[0]*(1e0/x-1e0/x_k)*(x_k)**2e0)) 
    for i in range(m):
        W = W + np.sum(np.where(dg[i+1]>0.,dg[i+1]*(x-x_k)*x_d[i],\
            -dg[i+1]*(1e0/x-1e0/x_k)*x_d[i]*x_k**2e0))
#
    return -W
#
# CONLIN: Dual function gradient
#
def dcon_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u)
#
    dW = np.zeros(m,dtype=np.float64)
#
    for i in range(m):
        dW[i] = dW[i] + g[i+1]
        dW[i] = dW[i] + np.sum(np.where(dg[i+1]>0.,dg[i+1]*(x-x_k),-dg[i+1]*(1./x-1./x_k)*(x_k)**2.))
#
#       for j in range(n):
#           if dg[i+1][j] > 0e0:
#               dW[i] = dW[i] + dg[i+1][j]*(x[j]-x_k[j])
#           else:
#               dW[i] = dW[i] - dg[i+1][j]*(1e0/x[j]-1e0/x_k[j])*(x_k[j])**2e0
#
    return -dW
#
