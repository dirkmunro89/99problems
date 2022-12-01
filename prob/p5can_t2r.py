#
import numpy as np
#
def init(g):
#
    n = 5; m = 1
    x_k = 5.*np.ones(n,dtype=float)
    x_l = 1.*np.ones(n,dtype=float)
    x_u = 10.*np.ones(n,dtype=float)
    x_d = np.ones(m,dtype=float)*1e6
    x_t='C'*n
#
    c_s=np.ones(m)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,x_t,x_d,c_s,aux
#
def apar(n):
#
    mov=0.2*np.ones(n)
    asf=[0.7,1.1]
#
    enf='none'
#
    kmx=100
    cnv=[1e-3,1e-3,1e-3,1e-3,1e-3]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_t, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[2e0*np.absolute(df_k[0])/x_k]
    for df in df_k[1:]: c_x.append((df[0],df[1],2e0*np.absolute(df[2])/x_k[df[1]]))
#
    L=x_k
    U=x_k
#
    d_l = np.maximum(x_k - mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k + mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    c1 = 0.0624
    c2 = np.array([61., 37., 19., 7., 1.], dtype=float)
#
    f = np.zeros((m + 1), dtype=float)
    f[0] =  c1*np.sum(x)
    f[1] = np.dot(c2, 1. / x ** 3.) - 1.
#
    df = [np.zeros(n, dtype=float)]
    df[0][:] = c1
#
    for i in range(n): df.append((0,i,-3.*c2[i]/x[i]**4.))
#
    return f, df
#
