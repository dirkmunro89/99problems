#
import osqp
import numpy as np
from scipy.optimize import minimize
from scipy import sparse
#
def t2dsqp(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x,c_s):
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
    A=sparse.csc_matrix(dg[1:])
    u=-g[1:]; l=-np.inf*np.ones(m,dtype=np.float64)
#   l=np.append(l,dx_l); u=np.append(u,dx_u)
#
    obj = osqp.OSQP()
    obj.setup(Q,J,A,l,u,verbose=True,warm_start=True,linsys_solver='mkl pardiso',\
        eps_abs=1e-6,eps_rel=1e-6,max_iter=int(1e6))
    res=obj.solve()
    sigma=1e1
    for k in range(10):
        a=np.where(   res.x-dx_u > 0  , 2.*(res.x-dx_u), 0)
        b=np.where(  -res.x+dx_l > 0  , -2.*(-res.x+dx_l), 0)
        c=np.where(  res.x-dx_u > 0  , 2, 0)
        d=np.where(  -res.x+dx_l > 0  , 2, 0)
        Jnew=J+(a+b)*sigma
        obj.update(q=Jnew)
        Q=sparse.csc_matrix((ddL+(c+d)*sigma, (ind, ind)), shape=(n, n))
        obj.update(Px=Q.data)
        res=obj.solve()
        sigma=sigma*10
    print(res.x-dx_u)
    print(-res.x+dx_l)
    print(res.x)
    stop
#
    if res.info.status != 'solved': print('WARNING: '+res.info.status)
#
    x_d[:]=res.y[:m]
    x=x_k+np.maximum(np.minimum(res.x,dx_u),dx_l)
#
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
#
    return x,x_d,q_k
#
