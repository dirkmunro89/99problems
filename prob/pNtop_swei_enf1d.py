#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
# specify subsolver here
#
from subs.mma87dual import mma87 as subs
#
# specify problem and algorithmic parameters here
#
def apar(n):
#   
    mov=1e-1*np.ones(n,dtype=float)
    asf=[0.7,1.2]
#
    enf='none'
#
    kmx=10
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, df, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.ones_like(df)*1e-6 
    c_x[1:]=0e0
#
    if k<=1:
        L = x_k-0.5*(x_u-x_l)
        U = x_k+0.5*(x_u-x_l)
    else:
        osc=(x_k-x_1)*(x_1-x_2)
        fac=np.ones_like(x_k)
        fac[np.where(osc>0)] = asf[1]
        fac[np.where(osc<0)] = asf[0]
        L=x_k-fac*(x_1-L_k)
        U=x_k+fac*(U_k-x_1)
#
    L_l = x_k - 10.*(x_u-x_l)
    L_u = x_k - 0.01*(x_u-x_l)
    U_l = x_k + 0.01*(x_u-x_l)
    U_u = x_k + 10.*(x_u-x_l)
    L=np.maximum(np.minimum(L, L_u), L_l)
    U=np.minimum(np.maximum(U, U_l), U_u)
#
    d_l = np.maximum(x_l, np.maximum(L + 0.1*(x_k-L), x_k - mov*(x_u-x_l)))
    d_u = np.minimum(x_u, np.minimum(U - 0.1*(U-x_k), x_k + mov*(x_u-x_l)))
#
    return c_x,mov,L,U,d_l,d_u
#
def init(g):
#
    mm=3
    nelx=20*mm
    nely=20*mm
    v_l = 0.1
    v_0 = 0.6
    v_u = 1.0
#
    ft = 1
    rmin = 1.1*mm
    dext=0
    felx = nelx+dext
    fely = nely+2*dext
#
    xPadd=np.zeros((felx,fely),dtype=float)
    tmp=-1*np.ones((nelx,nely),dtype=float)
    xPadd[:nelx,dext:fely-dext]=tmp
    xPadd[nelx-dext:,nely+dext:]=1.
    xPadd[nelx:,nely:nely+dext]=1.
    xPadd=xPadd.flatten()
#
    # BC's and support
    ndof=2*(nelx+1)*(nely+1)
    dofs=np.arange(2*(nelx+1)*(nely+1))
    fix=np.union1d(dofs[0:2*(nely+1):2],np.array([ndof-2,ndof-1]))
#
    # Set load
    frc=[(1,0)]
#
    pen = 3.0
    qen = 1.0
    muc = 1e-2
    Emin = 0e0; Emax=1.0
    gv = -9.81/nelx/nely
#
    n = nelx*nely
    m = 1
    x_l = np.ones(n,dtype=float)*1e-3
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,qen,muc,Emin,Emax,gv,g)
#
    return n,m,x_l,x_u,x_k,aux
#
def simu(n,m,x,aux,g):
#
    v_l=aux[2]
    v_u=aux[4]
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    [c,dc,v,dv]=topo2d_simu(n,m,x,aux,g)
#
    f[0]=c/360#0
#   f[1]=v/n/v_u-1.
    f[1]=-v/n/v_l+1.
#
    df[0][:] = dc/360#0
#   df[1][:] = dv/n/v_u
    df[1][:] = -dv/n/v_l
#
    return f, df
#
