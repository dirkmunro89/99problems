#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
def init(g):
#
    mm=6
    nelx=2*20*mm
    nely=20*mm
    v_l = 0.2
    v_0 = 0.2
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
    fix=np.union1d(dofs[0:2*(nely+1):2],np.array([ndof-1]))

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
