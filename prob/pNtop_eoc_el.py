#
import numpy as np
from prob.util.topo2d_el import topo2d_init
from prob.util.topo2d_el import topo2d_simu
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=2e0*np.absolute(df_k)/x_k
#
    L=x_k
    U=x_k
#
    d_l = np.maximum(x_k - mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k + mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def apar(n):
#
    mov=0.2*np.ones(n)
    asf=[0.7,1.1]
#
    enf='none'
#
    kmx=10
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def init(g):
#
    nelx = 180
    nely = 60
    v_l = 0.0
    v_0 = 0.4
    v_u = 0.4
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
    fix=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))

    # Set load
    frc=[(1,-1)]
#
    pen =  3.0
    qen =  0.0
    muc =  0.0
    Emin=1e-9; Emax=1.0
#
    gv=0.
    bf=np.array([[0.],[gv],[0.],[gv],[0.],[gv],[0.],[gv]])
    ei=np.array([[1.],[1.],[0.],[1.],[1.],[0.],[1.],[1.],[0.],[1.],[1.],[0.]])*0.
#
    n = nelx*nely
    m = 1
    x_l = np.zeros(n,dtype=float)
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,qen,muc,Emin,Emax,bf,ei,g)
#
    c_s=np.ones(m)
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def simu(n,m,x,aux,g):
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    [c,dc,v,dv]=topo2d_simu(n,m,x,aux,g)
#
    v_l=aux[2]
    v_u=aux[4]
#
    f[0]=c
    f[1]=v/n-v_u
#
    df[0][:] = dc
    df[1][:] = dv/n
#
    return f, df
#
