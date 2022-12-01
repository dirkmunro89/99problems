#
import cplex
import numpy as np
#
def t2c(n,m,x_k,x_t,x_d,d_l,d_u,g,dg,L,U,c_x,c_s):
#
    ddL = c_x[0].copy()
    for c in c_x[1:]: ddL[c[1]]=ddL[c[1]]+ x_d[c[0]]*c[2]
    ddL=np.maximum(np.absolute(ddL),1e-6)
#
    lb=d_l-x_k
    ub=d_u-x_k
#
    sub=cplex.Cplex()
    sub.set_problem_type(sub.problem_type.LP)
    sub.parameters.threads=8
    sub.parameters.read.datacheck=0
    sub.set_results_stream(None)
    sub.set_log_stream(None)
#
    rhs=[]; sen=[]
    for j in range(m):
        rhs.append((j,-g[j+1]))
        if c_s[j]: sen.append((j,'L'))
        else: sen.append((j,'E'))
#
    sub.variables.add(obj=dg[0],lb=lb,ub=ub,types=x_t)
    sub.linear_constraints.add(names=[None]*m)
    sub.linear_constraints.set_coefficients(dg[1:])
    sub.linear_constraints.set_rhs(rhs)
    sub.linear_constraints.set_senses(sen)
#   sub.objective.set_quadratic(ddL) #  HAS TO BE AFTER FIRST DEF OF OBJ !!!!
#
    sub.solve()
#
#   set up and solve the infeasibility problem, if needed
    if sub.solution.get_status()==103: 
        sub.feasopt(sub.feasopt.all_constraints())
#   else solve the fixed MILP problem, to get dual variables
#   else:
#       sub.set_problem_type(sub.problem_type.fixed_MILP)
#       sub.solve()
#       x_d[:]=-np.array(sub.solution.get_dual_values())

#   if infeasible now, try feasopt again and skip getting duals
#   if sub.solution.get_status()==3 or sub.solution.get_status()==103: 
#       sub.feasopt(sub.feasopt.all_constraints())
#   if feasible solution to MI problem was found, just continue
#   elif sub.solution.get_status()==120:
#       dummy=1
#   if feasible now, solve the fixed MILP to get duals
#   else:
#       sub.set_problem_type(sub.problem_type.fixed_MILP)
#       sub.solve()
#
#   print(sub.solution.get_status(), np.amax(np.absolute(x_d)))
#   stop
#
    x=np.minimum(np.maximum(x_k+sub.solution.get_values(),d_l),d_u)
#
    q_k = g.copy(); q_k[0]=q_k[0]+np.dot(dg[0],x-x_k)+np.dot(c_x[0]/2.,(x-x_k)**2.)
    for df in dg[1:]: q_k[df[0]+1]=q_k[df[0]+1]+df[2]*(x[df[1]]-x_k[df[1]])
    for c in c_x[1:]: q_k[c[0]+1]=q_k[c[0]+1]+c[2]/2.*(x[c[1]]-x_k[c[1]])**2.
#
    if sub.solution.get_status()==120: sol=0
    else: sol=1
#
    return x,x_d,q_k,sol
