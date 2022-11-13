#
import numpy as np
from scipy.optimize import minimize
#
def t2d(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x,c_s):
#
    bds=[[0e0*c_s[i]-1e6*(1-c_s[i]),1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(qpq_dual,x_d,args=(n,m,x_k,g,dg,d_l,d_u, c_x[0], c_x[1:]), \
        jac=dqpq_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    if sol.status != 0 or sol.success != True : print('Warning; subproblem')
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u, c_x[0], c_x[1:])
#
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
#
    return x,d,q_k
#
# QPQC: x in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    ddL = c0
    for c in cj:
        ddL[c[1]]=ddL[c[1]]+ x_d[c[0]]*c[2]
    ddL=np.maximum(np.absolute(ddL)),1e-6)
    tmp=dg[0]
    for df in dg[1:]
        tmp[df[1]]=tmp[df[1]]+x_d[df[0]]*df[2]
#
    return np.maximum(np.minimum(x_k - tmp/ddL, dx_u),dx_l)
#
# QPQC: Dual function value
#
def qpq_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj)
#
    ddL = c0
    for c in cj:
        ddL[c[1]]=ddL[c[1]]+ x_d[c[0]]*c[2]
    ddL=np.maximum(np.absolute(ddL)),1e-6)
#
    W=g[0]+np.dot(ddL/2.,(x-x_k)**2.)+np.dot(dg[0],x-x_k)
    for df in dg[1:]:
        W[df[1]]=W[df[1]]+x_d[df[0]]*(g[df[0]+1]+df[2]*(x[df[1]]-x_k[df[1]]))
#
    return -W
#
# QPQC: Dual gradient
#
def dqpq_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, c0, cj)
#
    dW=g[1:]+np.dot(cj/2e0,(x-x_k)**2e0) 
    for df in dg[1:]:
        dW[df[1]]=dW[df[1]]+df[2]*(x[df[1]]-x_k[df[1]]) 
#
    return -dW
#

