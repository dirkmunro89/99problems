#
import numpy as np
#
def init(g):
#
    n = 5; m = 1
    x_k = 5.*np.ones(n,dtype=float)
    x_l = 1.*np.ones(n,dtype=float)
    x_u = 1.e1*np.ones(n,dtype=float)
#
    c_s=np.ones(m)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#
    mov=0.2*np.ones(n)
    asf=[0.7,1.1]
#
    enf='none'
#
    kmx=14
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    if k<=1:
        L = x_k-mov*(x_u-x_l)
        U = x_k+mov*(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - asf[0]*(x_1 - L_k), x_k - asf[1]*(x_1 - L_k))
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + asf[0]*(U_k - x_1), x_k + asf[1]*(U_k - x_1))
#
    L = np.maximum(L, x_k-1e+2*(x_u-x_l))
    L = np.minimum(L, x_k-1e-5*(x_u-x_l))
    U = np.maximum(U, x_k+1e-5*(x_u-x_l))
    U = np.minimum(U, x_k+1e+2*(x_u-x_l))
#
    d_l = np.maximum(L+0.1*(x_k - L), x_l)
    d_u = np.minimum(U-0.1*(U - x_k), x_u)
#
    d_l= np.maximum(d_l, x_k-mov*(x_u-x_l))
    d_u= np.minimum(d_u, x_k+mov*(x_u-x_l))
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    c1 = 0.0624
    c2 = np.array([61, 37, 19, 7, 1], dtype=float)
#
    f = np.zeros((m + 1), dtype=float)
    f[0] = c1 * np.sum(x)
    f[1] = np.dot(c2, 1 / x ** 3) - 1
#
    df = np.zeros((m + 1, n), dtype=float)
    df[0][:] = c1
    df[1][:] = -3 * c2 / x ** 4
#
    return f, df
#


