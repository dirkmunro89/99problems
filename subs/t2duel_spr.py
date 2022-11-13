#
import numpy as np
from scipy.optimize import minimize
#
def t2d(n,m,x_k,x_d_k,d_l,d_u,g,dg,L,U,c_x,c_s):
#
    ddL = c_x[0]
    for c in c_x[1:]: ddL[c[1]]=ddL[c[1]]+ x_d_k[c[0]]*c[2]
    ddL=np.maximum(np.absolute(ddL),1e-6)
#
    bds=[[0e0*c_s[i]-1e6*(1-c_s[i]),1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(qp_dual,x_d_k.copy(),args=(n,m,x_k,g,dg,d_l,d_u,ddL), \
        jac=dqp_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    if sol.status != 0 or sol.success != True : print('Warning; subproblem')
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u, ddL)
#
    q_k = g; q_k[0]=q_k[0]+np.dot(dg[0],x-x_k)+np.dot(c_x[0]/2.,(x-x_k)**2.)
    for df in dg[1:]: q_k[df[0]]=q_k[df[0]]+df[2]*(x[df[1]]-x_k[df[1]])
    for c in c_x[1:]: q_k[c[0]]=q_k[c[0]]+c[2]/2.*(x[c[1]]-x_k[c[1]])**2.
#
    return x,d,q_k
#
# QP: x in terms of dual variables 
#
def x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    tmp=dg[0]
    for df in dg[1:]: tmp[df[1]]=tmp[df[1]]+x_d[df[0]]*df[2]
#
    return np.maximum(np.minimum(x_k - tmp/ddL, dx_u),dx_l)
#
# QP: Dual function value
#
def qp_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL)
#
    W=g[0]+np.dot(ddL/2.,(x-x_k)**2.)+np.dot(dg[0],x-x_k)
    for df in dg[1:]: W=W+x_d[df[0]]*(g[df[0]+1]+df[2]*(x[df[1]]-x_k[df[1]]))
#
    return -W
#
# QP: Dual gradient
#
def dqp_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL):
#
    x=x_dual(x_d, n, m, x_k, g, dg, dx_l, dx_u, ddL)
#
    dW=g[1:]
    for df in dg[1:]: dW[df[0]]=dW[df[0]]+df[2]*(x[df[1]]-x_k[df[1]]) 
#
    return -dW
#

