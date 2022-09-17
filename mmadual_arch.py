#
import numpy as np
from scipy.optimize import minimize
#
def mma(n,m,k,x_k,x_d,x_l,x_u,g,dg,x_1,x_2,m_r,m_a,L_k,U_k,a_f):
#
    if k<=1:
        L = x_k-m_a*(x_u-x_l)
        U = x_k+m_a*(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - a_f[0]*(x_1 - L_k), x_k - a_f[1]*(x_1 - L_k))
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + a_f[0]*(U_k - x_1), x_k + a_f[1]*(U_k - x_1))
#
    L = np.minimum(np.maximum(x_k-1e-2*(x_u-x_l),L),x_k-1e1*(x_u-x_l))
    U = np.minimum(np.maximum(x_k+1e-2*(x_u-x_l),U),x_k+1e1*(x_u-x_l))
#
    d_l = np.maximum(np.maximum(L + 0.1*(x_k-L), x_k - m_a*(x_u-x_l)),x_l)
    d_u = np.minimum(np.minimum(U - 0.1*(U-x_k), x_k + m_a*(x_u-x_l)),x_u)
#
    r = np.zeros((m+1),dtype=np.float64)
    p = np.zeros((m+1,n),dtype=np.float64)
    q = np.zeros((m+1,n),dtype=np.float64)
#
    for i in range(m+1):
        tmp=dg[i]
        p[i][:] = np.where(tmp>0e0, tmp*(U-x_k)**2e0, 0e0)
        q[i][:] = np.where(tmp>0e0, 0e0, -tmp*(x_k-L)**2e0)
        r[i] = g[i] - np.sum(p[i]/(U-x_k)) - np.sum(q[i]/(x_k-L))
#
    bds=[[1e-6,1e6] for i in range(m)]; tup_bds=tuple(bds)
    sol=minimize(mma_dual,x_d,args=(n,m,r,p,q,d_l,d_u,L,U), \
        jac=dmma_dual,method='L-BFGS-B',bounds=tup_bds, options={'disp':False})
#
    d=sol.x
    x=x_dual(d, n, m, r, p, q, d_l, d_u, L, U)
#
    return x,d,L,U
#
# primal variables in terms of dual variables 
#
def x_dual(x_d, n, m, r, p, q, d_l, d_u, L, U):
#
    tmp1= np.sqrt(np.maximum(p[0] + np.dot(x_d,p[1:]), np.zeros(n,dtype=np.float64)))
    tmp2= np.sqrt(np.maximum(q[0] + np.dot(x_d,q[1:]), np.zeros(n,dtype=np.float64)))
    x = np.minimum( np.maximum( (tmp1*L + tmp2*U) /  (tmp1 + tmp2), d_l), d_u )
#
    return x
#
# Dual function value
#
def mma_dual(x_d, n, m, r, p, q, d_l, d_u, L, U):
#
    x=x_dual(x_d, n, m, r, p, q, d_l, d_u, L, U)
    W=r[0]+np.sum(p[0]/(U-x)+np.dot(x_d,p[1:])/(U-x))+np.sum(q[0]/(U-x)+np.dot(x_d,q[1:])/(x-L)) \
        -np.dot(x_d,-r[1:])
#
    return -W
#
# Dual gradient
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
