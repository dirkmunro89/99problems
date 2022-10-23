#
from cvxopt import matrix, solvers, spdiag
import numpy as np
from scipy import sparse
#
def t2cvxo(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    dx_l=np.ones(n,dtype=np.float64)
    dx_u=np.ones(n,dtype=np.float64)
#
    ddL=np.zeros(n,dtype=np.float64)
#
    dx_l=d_l-x_k
    dx_u=d_u-x_k
#
    ddL=np.maximum(c_x[0]+np.dot(x_d.transpose(),c_x[1:]),0e0)
#
    J=dg[0]; ind = np.array(range(n))
    Q=sparse.csc_matrix((ddL, (ind, ind)), shape=(n, n))
    tmp=np.zeros((n,n),dtype=np.float64); np.fill_diagonal(tmp,1e0)
    A=np.append(np.append(dg[1:],tmp,axis=0),-tmp,axis=0)
    u=-g[1:]; l=-np.inf*np.ones(m,dtype=np.float64)
    l=np.append(l,dx_l); u=np.append(u,dx_u)
    h=-g[1:]
    h=np.append(h,dx_u)
    h=np.append(h,-dx_l)
#
    P = matrix(spdiag(list(ddL)),tc='d')
    q = matrix(J,tc='d')
    G = matrix(A,tc='d')
    h = matrix(h,tc='d')
#
    solvers.options['show_progress']=False
    sol=solvers.qp(P,q,G,h)
#
    x_d[:]=np.array(sol['s'][:m]).flatten()
    x=x_k+np.maximum(np.minimum(np.array(sol['x']).flatten(),dx_u),dx_l)
#
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
#
    return x,x_d,q_k
#
