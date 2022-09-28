#
import numpy as np
#
def init():
#
    n = 2; m = 2
    x_l = np.array([0.2, 0.1])
    x_u = np.array([4.0, 1.6])
    x_k = np.array([1.5, 0.5])
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,aux
#
def apar():
#   
    mov=2.
    asf=[0.7,1.1]
#
    enf='None'
#       
    kmx=8
    cnv=[1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, dg, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.ones_like(dg)*1e-6
#
    d_l = np.maximum(x_k/mov,x_l)
    d_u = np.minimum(mov*x_k,x_u)
#
    L=L_k
    U=U_k
#
    return c_x,L,U,d_l,d_u
#
def simu(n,m,x,aux):
#
    g = np.zeros((m + 1), dtype=float)
    dg = np.zeros((m + 1, n), dtype=float)
#
    c1 = 1.0; c2 = 0.124
    tmp1 = np.sqrt(1 + x[1] ** 2)
    tmp2 = 8 / x[0] + 1 / (x[0] * x[1])
    tmp3 = 8 / x[0] - 1 / (x[0] * x[1])
#
    g[0] = c1 * x[0] * tmp1
    g[1] = c2 * tmp1 * tmp2 - 1
    g[2] = c2 * tmp1 * tmp3 - 1
#
    tmp1 = np.sqrt(1 + x[1] ** 2)
    tmp2 = 8 / x[0] + 1 / (x[0] * x[1])
    tmp3 = 8 / x[0] - 1 / (x[0] * x[1])
    tmp4 = 2 * x[1]
#
    dg[0][0] = tmp1
    dg[0][1] = x[0] / (2 * tmp1) * tmp4
    dg[1][0] = -c2 * tmp1 * (8 / x[0] ** 2 + 1 / (x[0] ** 2 * x[1]))
    dg[1][1] = c2 / (2 * tmp1) * tmp4 * tmp2 - c2 * tmp1 / (x[0] * x[1] ** 2)
    dg[2][0] = -c2 * tmp1 * (8 / x[0] ** 2 - 1 / (x[0] ** 2 * x[1]))
    dg[2][1] = c2 / (2 * tmp1) * tmp4 * tmp3 + c2 * tmp1 / (x[0] * x[1] ** 2)
#
    return g, dg
#
