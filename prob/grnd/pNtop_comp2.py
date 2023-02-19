#
import numpy as np
from prob.util.ropo2d import topo2d_init
from prob.util.ropo2d import topo2d_simu
from cmls import t2rl as caml
from subs.t2dual import t2d as subs
#
def init(g):
#
    nelx = 180*1
    nely = 60*1
    v_l = 0.0
    v_0 = 0.25
    v_u = 0.25
#
    ft = 1
    rmin = 2
    dext=0#int(np.ceil(rmin))
    felx = nelx+dext
    fely = nely+2*dext
    xPadd = -1*np.ones(felx*fely,dtype=float)
#
    xPadd=np.zeros((felx,fely),dtype=float)
    tmp=-1*np.ones((nelx,nely),dtype=float)#.reshape(nelx,nely)
    xPadd[:nelx,dext:fely-dext]=tmp
    xPadd[nelx-dext:,nely+dext:]=1.
    xPadd[nelx:,nely:nely+dext]=1.
    xPadd=xPadd.flatten()
# 
    # BC's and support
    dofs=np.arange(2*(nelx+1)*(nely+1))
#   fix=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
    fix=[]#np.union1d(dofs[0:2*(nely+1):2])
    c=0
    for i in range(nelx+1):
        for j in range(nely+1):
            if j == 0:
                fix.append(c)
            c=c+1
            if j == 0:
                fix.append(c)
            c=c+1
    fix=np.array(fix)
    # Set load
    frc=[( int( (nely+1)*(nelx+1) + (nely+1) - 2), 1) , (2*(nely+1)-2,0.0), (2*(nely+1)*(nelx+1)-2,0.0)  ]
#
    pen =  3.0
    qen =  1.0
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
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,qen,muc,Emin,Emax,gv,g)
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
    f[1]=v/n/v_u-1.
#
    df[0][:] = dc
    df[1][:] = dv/n/v_u
#
    return f, df
#
