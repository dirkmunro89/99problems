#
import numpy as np
#
def init(g):
#
    n = 2; m = 2
    x_l = np.array([0.2, 0.1])
    x_u = np.array([4.0, 1.6])
    x_k = np.array([1.5, 0.5])
#
    c_s=np.ones(m)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#   
    mov=2.*np.ones(n)
    asf=[0.5,1.5]
#       
    enf='none'
#
    kmx=6
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    if k<=1:
        L = x_k-(x_u-x_l)
        U = x_k+(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - asf[0]*(x_1 - L_k), x_k - asf[1]*(x_1 - L_k))
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + asf[0]*(U_k - x_1), x_k + asf[1]*(U_k - x_1))
#
    i=0
    a_t=0.2
    L[i]=a_t*x_k[i]
    U[i]=x_k[i]/a_t
#
    d_l = np.maximum(np.maximum(x_k/mov, 1.01*L),x_l)
    d_u = np.minimum(np.minimum(mov*x_k, 0.99*U),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    c1 = 1.0; c2 = 0.124
    tmp1 = np.sqrt(1 + x[1] ** 2)
    tmp2 = 8 / x[0] + 1 / (x[0] * x[1])
    tmp3 = 8 / x[0] - 1 / (x[0] * x[1])
#
    f[0] = c1 * x[0] * tmp1
    f[1] = c2 * tmp1 * tmp2 - 1
    f[2] = c2 * tmp1 * tmp3 - 1
#
    tmp1 = np.sqrt(1 + x[1] ** 2)
    tmp2 = 8 / x[0] + 1 / (x[0] * x[1])
    tmp3 = 8 / x[0] - 1 / (x[0] * x[1])
    tmp4 = 2 * x[1]
#
    df[0][0] = tmp1
    df[0][1] = x[0] / (2 * tmp1) * tmp4
    df[1][0] = -c2 * tmp1 * (8 / x[0] ** 2 + 1 / (x[0] ** 2 * x[1]))
    df[1][1] = c2 / (2 * tmp1) * tmp4 * tmp2 - c2 * tmp1 / (x[0] * x[1] ** 2)
    df[2][0] = -c2 * tmp1 * (8 / x[0] ** 2 - 1 / (x[0] ** 2 * x[1]))
    df[2][1] = c2 / (2 * tmp1) * tmp4 * tmp3 + c2 * tmp1 / (x[0] * x[1] ** 2)
#
    return f, df
#
