#
import numpy as np
from subs.t2duel_spr import t2d as subs
#
def init(g):
#
    n = 5; m = 8
#   3 areas followed by 2 displacements
    x_l = np.array([0.0, 0.0, 0.0,-1.0,-1.0])
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
    kmx=1000
    cnv=[1e-3,1e-3,1e3,1e-3,1e-3]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, x_d, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[np.zeros_like(x_k)]
    if k > 0:
        sph = f_1 - f_k
        for df in df_k:
            sph[df[1]]=sph[df[1]]-df[2]*(x_1[df[1]]-x_k[df[1]])
        sph=sph/np.maximum(np.linalg.norm(x_1-x_k)**2./1e-6)
        print(sph)
        stop
#       for j in range(len(f_k)): c_x[j]=sph[j]
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
    df = [np.zeros(n,dtype=float)]
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
    df.append((0,3,-E*c[0]/l[0]**2./siy))
    df.append((0,4,E*100./l[0]**2./siy))
#
    df.append((1,3,-E*c[1]/l[1]**2./siy))
    df.append((1,4,E*100./l[1]**2./siy))
#
    df.append((2,3,-E*c[2]/l[2]**2./siy))
    df.append((2,4,E*100./l[2]**2./siy))
#
    df.append((3,3,E*c[0]/l[0]**2./siy))
    df.append((3,4,-E*100./l[0]**2./siy))
#
    df.append((4,3,E*c[1]/l[1]**2./siy))
    df.append((4,4,-E*100./l[1]**2./siy))
#
    df.append((5,3,E*c[2]/l[2]**2./siy))
    df.append((5,4,-E*100./l[2]**2./siy))
#
    df.append((6,0,E/l[0]**3*x[3] + E/l[0]**3.*x[4]))
    df.append((6,2,E/l[2]**3*x[3] - E/l[2]**3.*x[4]))
    df.append((6,3,k11/1e4))
    df.append((6,4,k12/1e4))
#
    df.append((7,0,E/l[0]**3.*x[3] + E/l[0]**3.*x[4]))
    df.append((7,1,E/l[1]**3*x[4]))
    df.append((7,2,-E/l[2]**3.*x[3] + E/l[2]**3.*x[4]))
    df.append((7,3,k12/1e4))
    df.append((7,4,k22/1e4))
#
    return f, df
#
