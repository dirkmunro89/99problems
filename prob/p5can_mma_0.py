#
import numpy as np
#
def init(g):
#
    n = 5; m = 1
    x_k = 5.*np.ones(n,dtype=float)
    x_l = 1.e-6*np.ones(n,dtype=float)
    x_u = 1.e6*np.ones(n,dtype=float)
#
    c_s=np.ones(m)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#
    mov=2.0*np.ones(n)
    asf=[0.5,1.50]
#
    enf='none'
#
    kmx=14
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    t=1e-6
#
    L=t*x_k
    U=x_k/t
#
    d_l = np.maximum(np.maximum(x_k/mov,1.01*L),x_l) #alpha
    d_u = np.minimum(np.minimum(mov*x_k,0.99*U),x_u) #beta
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


