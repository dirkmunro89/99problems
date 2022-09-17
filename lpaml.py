
from scipy.sparse import coo_matrix, csc_matrix
import numpy as np
import osqp

# linear program with move limit adaptation
def lpaml(n,x,volfrac_lo,volfrac_up,obj,dc,dv,xolder,xoldest,loop,mov_fct):
    prob = osqp.OSQP()
#   move limit adaptation
    if loop > 2:
        osc=(x-xolder)*(xolder-xoldest)/0.2
        for i in range(n):
            if osc[i] > 1e-9: mov_fct[i]=mov_fct[i]*1.2
            else: mov_fct[i]=mov_fct[i]*0.7
    dx_l=np.maximum(x-mov_fct*0.2,np.zeros(n))-x
    dx_u=np.minimum(x+mov_fct*0.2,np.ones(n))-x
    J=dc; tmp=np.zeros((n,n)); np.fill_diagonal(tmp,1e0)
    A=csc_matrix(np.append(dv.reshape((1,n)),tmp,axis=0))
    l=n*volfrac_lo-np.sum(x); u=n*volfrac_up-np.sum(x)
    l=np.append(l,dx_l); u=np.append(u,dx_u)
#   if loop==1:
    Q=csc_matrix((np.ones(n)*0e0,(np.array(range(n)),np.array(range(n)))),shape=(n,n))
    prob.setup(Q, J, A, l, u,verbose=False,warm_start=False)
#   else: prob.update(q=J, Ax=A.data, l=l, u=u)
    res=prob.solve()
    if res.info.status != 'solved': print('Solver warning: ', end='')
#   variable update
    xnew=x+np.maximum(np.minimum(res.x,dx_u),dx_l)
    return (xnew)

