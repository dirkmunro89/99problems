#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
# specify subsolver here
#
from subs.condual import con as subs
#
# specify problem and algorithmic parameters here
#
def apar(n):
#   
    mov=0.5*np.ones(n)
    asf=[0.7,1.1]
#  
    enf='none' # run with none and with t-r
#     
    kmx=200
    cnv=[1e-2,1e-2,1e-6,1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
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
#   run with and without L and U based aml
#
    d_l = np.maximum(np.maximum(x_k-mov*(x_u-x_l),x_l),L)
    d_u = np.minimum(np.minimum(x_k+mov*(x_u-x_l),x_u),U)
#   d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
#   d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def init(g):
#
    nelx=20
    nely=20
    v_l=0.01
    v_0=1.0
    v_u=0.8
#
    ft = 1
    rmin=1.1
    felx = nelx#+dext
    fely = nely#+2*dext
#
    xPadd=np.zeros((nelx,nely),dtype=float)
    tmp=-1*np.ones((nelx,nely),dtype=float)
    xPadd[:nelx,:fely]=tmp
    xPadd[nelx:,nely:]=1.
    xPadd[nelx:,nely:nely]=1.
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
    pen = 2.0
    muc = 0.
    Emin = 0e0; Emax=1.0
    gv = -9.81/800/1.
#
    n = nelx*nely
    m = 2
    x_l = np.ones(n,dtype=float)*2e-1
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,muc,Emin,Emax,gv,g)
#
    c_s=np.ones(m)
#
    return n,m,x_l,x_u,x_k,c_s,aux
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
    f[0]=c/n
    f[1]=v/n/v_u-1.
    f[2]=-v/n/v_l+1.
#
    df[0][:] = dc/n
    df[1][:] = dv/n/v_u
    df[2][:] = -dv/n/v_l
#
    return f, df
#
