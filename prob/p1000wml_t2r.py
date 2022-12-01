#
import numpy as np
from subs.t2cplx import t2c as subs
#
def init(g):
#
    n = 1000; m = 2
    x_k = 1e-5 * np.ones((n), dtype=float)
    x_l = 1e-6 * np.ones_like(x_k)
    x_u = 1e6 * np.ones_like(x_k)
    x_t = "C"*n
#
    aux=[]
#
    c_s=np.ones(m)
    x_d=np.ones(m,dtype=float)*1e6
#
    return n,m,x_l,x_u,x_k,x_t,x_d,c_s,aux
#
def apar(n):
#   
    mov=0.1*np.ones(n)
    asf=[0.7,1.1]
#  
    enf='none'
#     
    kmx=100
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_t, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[0e0*np.absolute(df_k[0])/x_k]
    for df in df_k[1:]: c_x.append((df[0],df[1],2e0*np.absolute(df[2])/x_k[df[1]]))
#
    L=np.zeros_like(x_k)
    U=np.zeros_like(x_k)
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    f=np.zeros((1+m),dtype=float)
    df = [np.zeros(n, dtype=float)]
#
    f[0]=np.sum(x)
    for i in range(950):
        f[1]=f[1]+1./x[i]
        f[2]=f[2]+1./x[i]
        df.append((0,i,-1./x[i]**2.))
        df.append((1,i,-1./x[i]**2.))
    for i in range(950,1000):
        f[1]=f[1]+(1.e-6)/x[i]
        f[2]=f[2]-(1.e-6)/x[i]
        df.append((0,i,-1e-6/x[i]**2.))
        df.append((1,i,1e-6/x[i]**2.))
#
    df[0]=np.ones_like(x,dtype=float)
#
    f[1]=f[1]-1000.
    f[2]=f[2]-900.
#
    return f, df
#
