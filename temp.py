#
import numpy as np
#
# specify subsolver here
#
from subs.t2dual import t2d as subs
#
# specify problem and algorithmic parameters here
#
def init():
#
    n = 1; m = 1
    x_l = np.array([1])
    x_u = np.array([1])
    x_k = np.array([1])
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,aux
#
def apar():
#   
    mov=1.
    asf=[0.7,1.1]
#       
    kmx=9
    cnv=[1e-6,1e-6]
#       
    return mov, asf, kmx, cnv
#
def caml(k, x_k, dg, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(dg)
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
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
    return g, dg
#
