#
import numpy as np
from subs.t2duel import t2d as subs
#
def init(g):
#
    n = 5; m = 8
#   3 areas followed by 2 displacements
#   x_l = np.array([0.2357, 0.1, 0.2357,-1.0,-1.0])
#   x_u = np.array([0.2357, 0.1, 0.2357, 1.0, 1.0])
#   x_k = np.array([0.2357, 0.1, 0.2357, 0., 0.])
    x_l = np.array([0.1, 0.1, 0.1,-1.0,-1.0])
    x_u = np.array([1., 1., 1., 1.0, 1.0])
    x_k = np.array([.5, .5, .5, 0., 0.])
#
    c_s = np.array([1,1,1,1,1,1,0,0]) # constraint sense, one is inequality
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#   
    mov=0.1*np.ones(n)
    asf=[0.7,1.1]
#
    enf='None'
#       
    kmx=100
    cnv=[1e-3,1e-3,1e3,1e-3,1e-3]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.ones_like(df_k)
    if k > 0:
        sph = 2.*(f_1 - f_k - np.dot(df_k,(x_1-x_k)))/np.maximum(np.linalg.norm(x_1-x_k)**2.,1e-6)
        for j in range(len(f_k)): c_x[j]=sph[j]
#
    d_l= np.maximum(x_l, x_k-mov*(x_u-x_l))
    d_u= np.minimum(x_u, x_k+mov*(x_u-x_l))
#
    L=L_k
    U=U_k
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    E = 30e6
    siy = 30e3
    c=np.array([100.,0.,-100.])
    l=np.sqrt(c**2. + 100.**2.)
    sig = E*(-x[3]*c+100.*x[4])/l**2.
    k11 = 100.**2.*E*(x[0]/l[0]**3. + x[2]/l[2]**3.)
    k12 = -100.**2.*E*(-x[0]/l[0]**3. + x[2]/l[2]**3.)
    k22 = 100.**2.*E*(x[0]/l[0]**3. + x[1]/l[1]**3. + x[2]/l[2]**3.)
#
    f[0] = 0.29 * np.sum(l*x[:3])/10e0
    f[1] = sig[0]/siy - 1.
    f[2] = sig[1]/siy - 1.
    f[3] = sig[2]/siy - 1.
    f[4] =-sig[0]/siy - 1.
    f[5] =-sig[1]/siy - 1.
    f[6] =-sig[2]/siy - 1.
    f[7] = k11*x[3]/10e3 + k12*x[4]/10e3 - 1.
    f[8] = k12*x[3]/10e3 + k22*x[4]/10e3
#
    df[0][0] = 0.29 * l[0]
    df[0][1] = 0.29 * l[1]
    df[0][2] = 0.29 * l[2]
    df[0]=df[0]/10.
#
#   df[1][0] = -E*x[3]/l[0]**2./siy
    df[1][3] = -E*c[0]/l[0]**2./siy
    df[1][4] = E*100./l[0]**2./siy
#
#   df[2][1] = -E*x[3]/l[1]**2./siy
    df[2][3] = -E*c[1]/l[1]**2./siy
    df[2][4] = E*100./l[1]**2./siy
#
#   df[3][2] = -E*x[3]/l[2]**2./siy
    df[3][3] = -E*c[2]/l[2]**2./siy
    df[3][4] = E*100./l[2]**2./siy
#
    df[4]=-df[1]
    df[5]=-df[2]
    df[6]=-df[3]
#
    df[7][0] = 100.**2.*E/l[0]**3*x[3] + 100.**2*E/l[0]**3.*x[4]
    df[7][1] = 0.
    df[7][2] = 100.**2.*E/l[2]**3*x[3] - 100.**2*E/l[2]**3.*x[4]
    df[7][3] = k11
    df[7][4] = k12
#
    df[8][0] = 100.**2*E/l[0]**3.*x[3] + 100.**2*E/l[0]**3.*x[4]
    df[8][1] = 100.**2.*E/l[1]**3*x[4]
    df[8][2] = -100.**2.*E/l[2]**3.*x[3] + 100.**2.*E/l[2]**3.*x[4]
    df[8][3] = k12
    df[8][4] = k22
#
    df[7]=df[7]/10e3
    df[8]=df[8]/10e3
#
    return f, df
#
