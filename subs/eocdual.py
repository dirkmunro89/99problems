#
import numpy as np
from scipy.optimize import minimize
#
def eoc(n,m,k,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    if m > 1: print('ERROR'); stop
#
    bds=[[1e-6,1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(eoc_dual,x_d,args=(n,m,x_k,g,dg,d_l,d_u), \
        jac=deoc_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    d=sol.x
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u)
#
    return x,d
#
# primal variables in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, d_l, d_u):
#
    beta=-(dg[0]*x_k**2e0)/x_d[0]/dg[1]
    x = np.maximum(np.minimum(beta**(1e0/2e0),d_u),d_l)
#
    return x
#
# Dual function value
#
def eoc_dual(x_d, n, m, x_k, g, dg, d_l, d_u):
#
    x=x_dual(x_d, n, m, x_k, g, dg, d_l, d_u)
#
    W = g[0] + x_d[0]*g[1] + x_d[0]*np.dot(dg[1],x-x_k)
    W = W + np.dot((x - x_k)*(x_k/(x+1e-6)),dg[0])
#
    return -W
#
# Dual gradient
#
def deoc_dual(x_d, n, m, x_k, g, dg, d_l, d_u):
#
    x=x_dual(x_d, n, m, x_k, g, dg, d_l, d_u)
#
    dW=g[1]+np.dot(dg[1],(x-x_k))
#
    return -dW
#
