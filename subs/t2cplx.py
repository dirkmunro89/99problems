#
import time
import cplex
import numpy as np
from scipy.optimize import minimize
from scipy import sparse
from docplex.mp.sktrans.transformers import CplexTransformer
#
def t2cplx(n,m,x_k,x_d,d_l,d_u,g,dg,L,U,c_x):
#
    dx_l=np.ones(n,dtype=np.float64)
    dx_u=np.ones(n,dtype=np.float64)
#
    ddL=np.zeros(n,dtype=np.float64)
#
    dx_l=d_l-x_k
    dx_u=d_u-x_k
#
    ddL=np.maximum(c_x[0]+np.dot(x_d.transpose(),c_x[1:]),1e-6)
#
    prob=cplex.Cplex()
    prob.set_results_stream(None)
    prob.set_log_stream(None)
    prob.parameters.threads=8
    prob.parameters.read.datacheck=0
    prob.variables.add(obj=dg[0][:],lb=dx_l,ub=dx_u)
    ind = range(n)#[i for i in range(n)]
    lin_expr=[]
    for j in range(m): lin_expr.append( cplex.SparsePair( ind = ind, val=dg[j+1] ))
    rhs=-g[1:]
    t1=time.time()
    prob.linear_constraints.add(lin_expr=lin_expr,rhs=rhs, senses=['L']*m)#['L' for i in range(m)])
    t2=time.time()
    print(t2-t1)
    prob.objective.set_quadratic(ddL)
    prob.solve()
    if prob.solution.get_status() == 3:
        prob.feasopt(prob.feasopt.linear_constraints())
    else:
        x_d[:]=-np.array(prob.solution.get_dual_values())
#
    x=x_k+np.maximum(np.minimum(prob.solution.get_values(),dx_u),dx_l)
#
    q_k = g+np.dot(dg,x-x_k)+np.dot(c_x/2.,(x-x_k)**2.)
#
    return x,x_d,q_k
#
