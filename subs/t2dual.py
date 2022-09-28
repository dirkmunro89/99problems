#
import numpy as np
from scipy.optimize import minimize
#
def t2d(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    bds=[[1e-6,1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(qpq_dual,x_d,args=(n,m,x_k,g,dg,d_l,d_u, c_x[0], c_x[1:]), \
        jac=dqpq_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    if sol.status != 0 or sol.success != True : print('Warning; subproblem')
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u, c_x[0], c_x[1:])
#
#   tmp=-np.dot(dg,x-x_k)-np.dot(c_x/2.,(x-x_k)**2.)
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
    dq=g[0]-q_k[0]
#   print(tmp,dq)
#
    return x,d,dq,q_k
#
# QPQC: x in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    ddL=(c0 + np.dot(x_d,cj))
    tmp=(dg[0]+np.dot(x_d,dg[1:]))
#
    return np.maximum(np.minimum(x_k - tmp/ddL, dx_u),dx_l)
#
# QPQC: Dual function value
#
def qpq_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj)
#
    ddL=(c0 + np.dot(x_d,cj))
#
    W=g[0]+np.dot(dg[0],x-x_k)+np.dot(ddL/2.,(x-x_k)**2.)+np.dot(x_d,(g[1:]+np.dot(dg[1:],(x-x_k))))
#
    return -W
#
# QPQC: Dual gradient
#
def dqpq_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj)
#
    dW=g[1:]+np.dot(dg[1:],(x-x_k)) + np.dot(cj/2e0,(x-x_k)**2e0)
#
    return -dW
#

