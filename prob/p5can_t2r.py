#
import numpy as np
#
def init(g):
#
    n = 5; m = 1
    x_k = 5.*np.ones(n,dtype=float)
    x_l = 1.*np.ones(n,dtype=float)
    x_u = 10.*np.ones(n,dtype=float)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,aux
#
def apar():
#
    mov=0.2
    asf=[0.7,1.1]
#
    enf='none'
#
    kmx=14
    cnv=[1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, df, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=2e0*np.absolute(df)/x_k
#
    L=x_k
    U=x_k
#
    d_l = np.maximum(x_k - mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k + mov*(x_u-x_l),x_u)
#
    return c_x,L,U,d_l,d_u
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


