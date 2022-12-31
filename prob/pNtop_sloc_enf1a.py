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
    mov=0.1e0*np.ones(n,dtype=float)
    asf=[0.5,1.5]
#
    enf='none'
#
    kmx=1000
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=2e0*np.absolute(df_k)/x_k
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
    mm=1
    nelx=40*mm
    nely=20*mm
    v_l = 0.01
    v_0 = 0.01
    v_u = 0.5
#
    ft = 1
    rmin = 2.*mm
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
    fix=np.union1d(dofs[0:2*(nely+1):2],np.array([ndof-1]))
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
#   for i in range(nelx):
#           c=i*nely+0
#           x_l[c]=1.
#           x_k[c]=1.
#           c=i*nely+1
#           x_l[c]=1.
#           x_k[c]=1.
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
    [nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,pen,qen,muc,Emin,Emax,gv,\
        ndof,KE,H,Hs,iK,jK,edofMat,fixed,free,f,u,fig,im]=aux
#
    f = np.zeros((m + 1), dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    [c,dc,v,dv]=topo2d_simu(n,m,x,aux,g)
#
    f[0]=c/n/36#np.sum(x)/n
    f[1]=v/n/v_u-1.
    f[2]=-v/n/v_l+1.
#
    xp=x#np.asarray(H*x[np.newaxis].T/Hs)[:,0]#[pad] ###
    c=0
    a=3.
    p=3.#/1.
    q=2.#/1.#/1.#1.1#p
    z=1.
    for i in range(1,felx-1):
        for j in range(fely-1):
            c=i*fely+j
            l=(i-1)*fely+j+1
            r=(i+1)*fely+j+1
            u=i*fely+j+1
            tmp=(xp[l]**(a*p)+xp[u]**(a*p)+xp[r]**(a*p))/z
            f[3+c]=xp[c]**q - tmp**(1./a) #- 1e-3#-10e0
            df[3+c][c]=q*xp[c]**(q-1)
            df[3+c][l]=  -(1./a)/z*(tmp)**(1./a-1.)*a*p*xp[l]**(a*p-1)
            df[3+c][r]=  -(1./a)/z*(tmp)**(1./a-1.)*a*p*xp[r]**(a*p-1)
            df[3+c][u]=  -(1./a)/z*(tmp)**(1./a-1.)*a*p*xp[u]**(a*p-1)
#
    df[0][:] = dc/n/36#np.ones(n,dtype=float)/n#dc/360#0
    df[1][:] = dv/n/v_u
    df[2][:] = -dv/n/v_l
#
#   for j in range(3,m+1):
#       df[j][:] = np.asarray(H*(df[j][np.newaxis].T/Hs))[:,0]
#
    return f, df
#
