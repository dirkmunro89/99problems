#
import numpy as np
from scipy.optimize import minimize
#
def dlp(n,m,k,x_k,x_d,x_l,x_u,g,dg,x_1,x_2,m_r,m_a,L_k,U_k,a_f):
#
    if k<=1:
        L = x_k-m_a*(x_u-x_l)
        U = x_k+m_a*(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - a_f[0]*(x_1 - L_k), x_k - (x_1 - L_k)*a_f[1])
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + a_f[0]*(U_k - x_1), x_k + (U_k - x_1)*a_f[1])
#
    L = np.minimum(np.maximum(x_k-1e-2*(x_u-x_l),L),x_k-1e1*(x_u-x_l))
    U = np.minimum(np.maximum(x_k+1e-2*(x_u-x_l),U),x_k+1e1*(x_u-x_l))
#
    d_l = np.maximum(np.maximum(L + 0.1*(x_k-L), x_k - m_a*(x_u-x_l)),x_l)
    d_u = np.minimum(np.minimum(U - 0.1*(U-x_k), x_k + m_a*(x_u-x_l)),x_u)
#
    ddg=np.zeros((m+1,n),dtype=float)
    ddg[0]=1e-6*np.ones(n,dtype=float) # try one
#
    bds=[[1e-6,1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(qpq_dual,x_d,args=(n,m,x_k,g,dg,d_l,d_u, ddg[0], ddg[1:]), \
        jac=dqpq_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    if sol.status != 0 or sol.success != True : print('Warning; subproblem')
#
    d=sol.x
#
    x=x_dual(d, n, m, x_k, g, dg, d_l, d_u, ddg[0], ddg[1:])
#
    return x,d,L,U
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

