#
#   todo
#
import numpy as np
from subs.t2osqp import t2osqp as subs
#from subs.t2dual import t2d as subs
#
def init(g):
#
    n = 1000; m = 1000
    x_k = 1e-5 * np.ones((n), dtype=float)
    x_l = 1e-6 * np.ones_like(x_k)
    x_u = 1e6 * np.ones_like(x_k)
#
    aux=[]
#
    c_s=np.ones(m)
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#   
    mov=0.1*np.ones(n)
    asf=[0.7,1.1]
#  
    enf='none'
#     
    kmx=1000
    cnv=[1e-5,1e-5,1e-12,1e-12,1e-12]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    L=np.zeros_like(x_k)
    U=np.zeros_like(x_k)
#
#   T2R
    c_x=-2e0*df_k/x_k
#
    c_x[0]=1e-6*np.ones_like(x_k)
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    f=np.zeros((1+m),dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    f=func(n,m,x)
#
    dx=1e-4
    for i in range(n):
        x0 = x[i]
        x[i] += dx
        fd = func(n,m,x)
        x[i] = x0
        df[:, i] = (fd - f) / dx
#
    np.savetxt('xxx.dat',x)
#
    return f, df
#
def func(n,m,x):
#
    f = np.zeros((m + 1), dtype=float)
#
    f[0]=np.sum(x)
    for i in range(949):
        f[1]=f[1]+1./x[i]
        f[2]=f[2]+1./x[i]
    for i in range(950,1000):
        f[1]=f[1]+(1e-6)/x[i]
        f[2]=f[2]-(1e-6)/x[i]
#
    f[1]=f[1]-1000.
    f[2]=f[2]-900.
#
    return f
#
