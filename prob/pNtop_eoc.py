#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
def caml(k, x_k, dg, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=0e0*np.absolute(dg)/x_k
#
    L=x_k
    U=x_k
#
    d_l = np.maximum(x_k - mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k + mov*(x_u-x_l),x_u)
#
    return c_x,L,U,d_l,d_u
#
def apar():
#
    mov=0.2
    asf=[0.7,1.1]
#
    enf='none'
#
    kmx=10
    cnv=[1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def init():
#
    nelx = 180
    nely = 60
    v_l = 0.0
    v_0 = 0.4
    v_u = 0.4
#
    vis=False
#
    ft = 1
    rmin = 5.4
    felx = 180
    fely = 60
    xPadd = -1*np.ones(felx*fely,dtype=float)
#
#   tmp=x.reshape(nelx,nely)
#   x_ext=np.zeros((nelx_ext,nely_ext))
#   x_ext[:nelx,d_ext:nely_ext-d_ext]=tmp
#   x_ext[nelx-d_ext:,nely+d_ext:]=1.
#   x_ext[nelx:,nely:nely+d_ext]=1.
#   x_ext=x_ext.flatten()
# 
    # BC's and support
    dofs=np.arange(2*(nelx+1)*(nely+1))
    fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))

    # Set load
    force=[(1,-1)]
#
    pen =  3.0
    muc =  0.0
    Emin=1e-9; Emax=1.0
    gv=0e0#-9.81#/nelx/nely
#
    n = nelx*nely
    m = 1
    x_l = np.zeros(n,dtype=float)
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fixed,force,pen,muc,Emin,Emax,gv,vis)
#
    return n,m,x_l,x_u,x_k,aux
#
def simu(n,m,x,aux):
#
    g = np.zeros((m + 1), dtype=float)
    dg = np.zeros((m + 1, n), dtype=float)
#
    [c,dc,v,dv]=topo2d_simu(n,m,x,aux)
#
    v_l=aux[2]
    v_u=aux[4]
#
    g[0]=c
    g[1]=v/n-v_u
#
    dg[0][:] = dc
    dg[1][:] = dv/n
#
    return g, dg
#
