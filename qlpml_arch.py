#
import numpy as np
from scipy.optimize import minimize
from scipy.sparse import coo_matrix, csc_matrix
import osqp
#
def qlp(n,m,k,x_k,x_d,x_l,x_u,g,dg,x_1,x_2,m_r,m_a,L_k,U_k,a_f):
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
    d_l = np.maximum(np.maximum(L + 0.1*(x_k-L), x_k - m_a*(x_u-x_l)),x_l) - x_k
    d_u = np.minimum(np.minimum(U - 0.1*(U-x_k), x_k + m_a*(x_u-x_l)),x_u) - x_k
#
    J=dg[0]; tmp=np.zeros((n,n)); np.fill_diagonal(tmp,1e0)
    tmp2=np.append(dg[2].reshape((1,n)),tmp,axis=0)
    tmp3=np.append(dg[1].reshape((1,n)),tmp2,axis=0)
#
    A=csc_matrix(tmp3)
#
    l=-np.ones(m)*1e12; u=-g[1:]#n*volfrac_up-np.sum(x)
    l=np.append(l,d_l); u=np.append(u,d_u)
#
    prob = osqp.OSQP()
    Q=csc_matrix((np.ones(n)*1e-6,(np.array(range(n)),np.array(range(n)))),shape=(n,n))
    prob.setup(Q, J, A, l, u,verbose=False,warm_start=False,eps_abs=1e-5,eps_rel=1e-5, \
        scaling=100,max_iter=10000)
#   else: prob.update(q=J, Ax=A.data, l=l, u=u)
    res=prob.solve()
    if res.info.status != 'solved': print('Solver warning: ', end='')
#   variable update
    x=x_k+np.maximum(np.minimum(res.x,d_u),d_l)
#
    return x,x_d,L,U
#
