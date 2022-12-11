#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
# specify subsolver here
#
from subs.t2dual import t2d as subs
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
    kmx=5000
    cnv=[1e-3,1e0,1e0,1e0,1e0]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    c_x[0]=1e-6
#
    if k>2:
        mov=np.where((x_k-x_1)*(x_1-x_2) <= 0., mov*asf[0], mov*asf[1])
#
    mov=np.minimum(np.maximum(mov,0.5e-3),.1)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def init(g):
#
    mm=3
    nelx=20*mm
    nely=20*mm
    v_l = 0.1
    v_0 = 1.0
    v_u = 1.0
#
    ft = 1
    rmin = 1.1*mm
    dext=0#int(np.ceil(rmin))
    felx = nelx+dext
    fely = nely+2*dext
#
    xPadd=np.zeros((felx,fely),dtype=float)
    tmp=-1*np.ones((nelx,nely),dtype=float)#.reshape(nelx,nely)
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
    muc = 0.0
    Emin = 1e-9; Emax=1.0
    gv = -1.
#
    n = nelx*nely
    m = 1
    x_l = np.ones(n,dtype=float)*0e0
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,qen,muc,Emin,Emax,gv,g)
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
    f[0]=c
    f[1]=-v/n/v_l+1.
#
    df[0][:] = dc
    df[1][:] = -dv/n/v_l
#
    return f, df
#
