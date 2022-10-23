#
import numpy as np
from prob.util.topo2d import topo2d_init
from prob.util.topo2d import topo2d_simu
#
# specify subsolver here
#
from subs.t2dual import t2d as subs
from subs.t2cvxo import t2cvxo as subs
from subs.t2osqp import t2osqp as subs
from subs.t2cplx import t2cplx as subs
#
# specify problem and algorithmic parameters here
#
def apar(n):
#   
    mov=0.5e0*np.ones(n,dtype=float)
    asf=[0.5,1.5]
#
    enf='t-r'
#
    kmx=100
    cnv=[1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, df, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=2e0*np.absolute(df)/x_k
#   c_x[1:]=0e0
#
    c_x[0]=np.maximum(c_x[0],1e-6)
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
    mm=2
    nelx=20*mm
    nely=20*mm
    v_l = 0.01
    v_0 = 1.0
    v_u = 1.0
#
    ft = 1
    rmin = 1.5#*mm
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
#   frc=[(1,0)]
    frc=[(1,-1)]
#
    pen = 3.0
    qen = 1.0
    muc = 1e-2
    Emin = 0e0; Emax=1.0
    gv = 0.#-9.81/nelx/nely
#
    n = nelx*nely
    m = 2 + n
    x_l = np.ones(n,dtype=float)*1e-6
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
    for i in range(nelx):
            c=i*nelx+0
            x_l[c]=1.
            x_k[c]=1.
            c=i*nelx+1
            x_l[c]=1.
            x_k[c]=1.
    aux=topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fix,frc,pen,qen,muc,Emin,Emax,gv,g)
#
    return n,m,x_l,x_u,x_k,aux
#
def simu(n,m,x,aux,g):
#
    v_l=aux[2]
    v_u=aux[4]
    [nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,pen,qen,muc,Emin,Emax,gv,\
        ndof,KE,H,Hs,iK,jK,edofMat,fixed,free,f,u,fig,im]=aux
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    [c,dc,v,dv]=topo2d_simu(n,m,x,aux,g)
#
    f[0]=np.sum(x)/n
    f[1]=np.sum(x)/n/v_u-1.
    f[2]=-np.sum(x)/n/v_l+1.
#
    c=0
    a=3.
    p=3.
    q=1./p
    for i in range(1,felx-1):
        for j in range(fely-1):
            c=i*felx+j
            l=(i-1)*felx+j+1
            r=(i+1)*felx+j+1
            u=i*felx+j+1
            tmp=(x[l]**(a*p)+x[u]**(a*p)+x[r]**(a*p))
            f[3+c]=x[c]**q - tmp**(1./a)
            df[3+c][c]=q*x[c]**(q-1)
            df[3+c][l]=  -(1./a)*tmp**(1./a-1.)*a*p*x[l]**(a*p-1)
            df[3+c][r]=  -(1./a)*tmp**(1./a-1.)*a*p*x[r]**(a*p-1)
            df[3+c][u]=  -(1./a)*tmp**(1./a-1.)*a*p*x[u]**(a*p-1)
#
    df[0][:] = np.ones(n,dtype=float)/n#dc/360#0
    df[1][:] = np.ones(n,dtype=float)/n/v_u
    df[2][:] = -np.ones(n,dtype=float)/n/v_l
#
    return f, df
#
