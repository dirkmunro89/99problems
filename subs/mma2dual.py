#
import numpy as np
from scipy.optimize import minimize
#
def mma(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    r = np.zeros((m+1),dtype=float)
    p = np.zeros((m+1,n),dtype=float)
    q = np.zeros((m+1,n),dtype=float)
#
    for i in range(m+1):
        tmp=dg[i]
        dpos=np.maximum(dg[i],0.)
        dneg=np.maximum(-dg[i],0.)
        p[i][:] = (1.001*dpos+0.001*dneg + 1e-5/(1.-1e-3))*(U-x_k)**2e0
        q[i][:] = (0.001*dpos+1.001*dneg + 1e-5/(1.-1e-3))*(x_k-L)**2e0
        r[i] = g[i] - np.sum(p[i]/(U-x_k) + q[i]/(x_k-L))
#
    bds=[[1e-6,1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(mma_dual,x_d,args=(n,m,r,p,q,d_l,d_u,L,U), \
        jac=dmma_dual,method='L-BFGS-B',bounds=tup_bds,options={'disp':False})
#
    d=sol.x
#
    x=x_dual(d, n, m, r, p, q, d_l, d_u, L, U)
#
    q_k=r+np.sum(p/(U-x),axis=1)+np.sum(q/(x-L),axis=1)
#
    return x,d,q_k
#
# primal variables in terms of dual variables 
#
def x_dual(x_d, n, m, r, p, q, d_l, d_u, L, U):
#
    tmp1= np.sqrt(np.maximum(p[0] + np.dot(x_d,p[1:]), np.zeros(n,dtype=float)))
    tmp2= np.sqrt(np.maximum(q[0] + np.dot(x_d,q[1:]), np.zeros(n,dtype=float)))
    x = np.minimum( np.maximum( (tmp1*L + tmp2*U) /  (tmp1 + tmp2), d_l), d_u )
#
    return x
#
# Dual function value
#
def mma_dual(x_d, n, m, r, p, q, d_l, d_u, L, U):
#
    x=x_dual(x_d, n, m, r, p, q, d_l, d_u, L, U)
    W=r[0]+np.sum(p[0]/(U-x)+np.dot(x_d,p[1:])/(U-x))+np.sum(q[0]/(x-L)+np.dot(x_d,q[1:])/(x-L)) \
        -np.dot(x_d,-r[1:])
#
    return -W
#
# Dual function gradient
#
def dmma_dual(x_d, n, m, r, p, q, d_l, d_u, L, U):
#
    x=x_dual(x_d, n, m, r, p, q, d_l, d_u, L, U)
    tmp11=np.sum(p[1:]/(U-x),axis=1)
    tmp22=np.sum(q[1:]/(x-L),axis=1)
    dW = r[1:] + tmp11 + tmp22
#
    return -dW
#
