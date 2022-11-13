#
import numpy as np
from scipy.optimize import minimize
#
def t2d(n,m,x_k,x_d_k,d_l,d_u,g,dg,L,U,c_x,c_s):
#
    ddL=np.maximum(np.absolute(c_x[0] + np.dot(x_d_k,c_x[1:])),1e-6)
#
    bds=[[0e0*c_s[i]-1e6*(1-c_s[i]),1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(qp_dual,x_d_k.copy(),args=(n,m,x_k,g,dg,d_l,d_u,ddL), \
        jac=dqp_dual,method='L-BFGS-B',bounds=tup_bds, \
        options={'disp':False,'gtol':1e-16,'ftol':1e-16,'maxls':100})
#
    if sol.status != 0 or sol.success != True : print('Warning; subproblem')
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u, ddL)
#
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
#
    return x,d,q_k
#
# QP: x in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    tmp=(dg[0]+np.dot(x_d,dg[1:]))
#
    return np.maximum(np.minimum(x_k - tmp/ddL, dx_u),dx_l)
#
# QP: Dual function value
#
def qp_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL)
#
    W=g[0]+np.dot(dg[0],x-x_k)+np.dot(ddL/2.,(x-x_k)**2.)+np.dot(x_d,(g[1:]+np.dot(dg[1:],(x-x_k))))
#
    return -W
#
# QP: Dual gradient
#
def dqp_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL)
#
    dW=g[1:]+np.dot(dg[1:],(x-x_k)) 
#
    return -dW
#

