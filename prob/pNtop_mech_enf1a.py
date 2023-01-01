#
import numpy as np
from prob.util.topo2d_sech import topo2d_init
from prob.util.topo2d_sech import topo2d_simu
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
    asf=[0.5,1.5]
#
    enf='none'
#
    kmx=1000
    cnv=[1e-1,1e-1,1e-3,1e-3,1e-3]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=2e0*np.absolute(df_k)/np.maximum(x_k,1e-6)
    c_x[1:]=0e0
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    if k>1e8:
        mov=np.where((x_k-x_1)*(x_1-x_2) <= 0., mov*asf[0], mov*asf[1])
        mov=np.minimum(np.maximum(mov,1e-3),0.1)
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
    mm=4
    nelx=20*mm
    nely=10*mm
    v_l = 0.1
    v_0 = 0.3
    v_u = 0.3
#
    ft = 1
    rmin = 2./4.*mm
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
    fix=np.union1d(dofs[1:2*(nely+1)*(nelx+1)-1:2*(nely+1)], dofs[2*(nely+1)-1:2*(nely+1)-5:-1])
#
    pen = 3.0
    qen = 1.0
    muc = 0e0
    Emin = 1e-9; Emax=1.0
    gv = 0.0/nelx/nely
#
    n = nelx*nely
    m = 1
    x_l = np.zeros(n,dtype=float)
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,pen,qen,muc,Emin,Emax,gv,g)
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
    f[0]=c#/360#0
    f[1]=v/n/v_u-1.
#   f[2]=-v/n/v_l+1.
#
    df[0][:] = dc#/360#0
    df[1][:] = dv/n/v_u
#   df[2][:] = -dv/n/v_l
#
    return f, df
#
