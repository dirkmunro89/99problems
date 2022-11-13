#
#   todo
#
import numpy as np
from subs.t2dual import t2d as subs
#
def init(g):
#
    n = 1000; m = 2
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
    kmx=100
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
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
    c_x=np.absolute(2e0*df_k/x_k)
#
    c_x[0]=0e0*np.ones_like(x_k)
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
    f[0]=np.sum(x)#/1e3
    for i in range(950):
        f[1]=f[1]+1./x[i]
        f[2]=f[2]+1./x[i]
        df[1][i]= -1./x[i]**2.
        df[2][i]= -1./x[i]**2.
    for i in range(950,1000):
        f[1]=f[1]+(1.e-6)/x[i]
        f[2]=f[2]-(1.e-6)/x[i]
        df[1][i]= -1.e-6/x[i]**2.
        df[2][i]=1.e-6/x[i]**2.
#
    df[0]=np.ones_like(x,dtype=float)#/1e3
#
#   df[1]=df[1]/1000.
#   df[2]=df[2]/900.
    f[1]=f[1]-1000.
    f[2]=f[2]-900.
#
    return f, df
#
