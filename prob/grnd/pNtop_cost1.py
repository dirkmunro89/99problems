#
import numpy as np
from prob.util.ropo2d import topo2d_init
from prob.util.ropo2d import topo2d_simu
#
from prob.util.cost2d import cost2d_simu
from prob.util.cost2d import plot
#
def init(g):
#
    nelx = 90#0#0*1
    nely = 30#30#0*1
    v_l = 0.0
    v_0 = 1.
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
    n = nelx*nely + 2
    m = 2
    x_l = np.zeros(n,dtype=float)+1e-3
    x_u = np.ones(n,dtype=float)
    x_k = v_0*np.ones(n,dtype=float)
#
    x_l[-2] = 0.
    x_l[-1] = 1.
    x_u[-2] = 0.
    x_u[-1] = 1.
    x_k[-2] = 0.
    x_k[-1] = 1.
#
#   x_l[-2] = -1.
#   x_l[-1] = -1.
#   x_u[-2] = 1.
#   x_u[-1] = 1.
#   x_k[-2] = 0.1
#   x_k[-1] = 1.
#
#
#   matrices for spatial gradients
#
#   [Gx,Gy] = grad_init(nelx,nely,2.) # put into aux
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
    [nelx,nely,_,_,_,_,rmin,_,_,_,_,_,_,_,_,_,\
        _,_,H,Hs,Gx,Gy,_,_,_,_,_,_,_,_]=aux
#
    x_k=np.zeros(nelx*nely,dtype=float)
    x_k[:]=x[:nelx*nely]
    b_k=x[-2:]
#
    [c,dc,v,dv]=topo2d_simu(nelx*nely,m,x_k,aux,g)
#
    v_l=aux[2]
    v_u=aux[4]
#
    f[0]=v/n
    f[1]=c/7.289e0/1.-1.
#
    df[0][:nelx*nely] = dv/n
    df[1][:nelx*nely] = dc/7.289e0/1.
#
#   cost
#
    ovrh = np.cos(np.pi/180.*50.)
#
    [srf, dsrfdb, dsrfdx, vol, dvoldx]=cost2d_simu(x_k, b_k, Gx, Gy, H, Hs)
#
#   f[1] = c/7.289e+00 - 1.
#   df[1][:nelx*nely] = dc/7.289e+00
#
    f[0] = srf#/143.
    df[0][:nelx*nely] = dsrfdx#/143.
    df[0][-2:] = dsrfdb#/143.
#
    f[2]=-v/n/0.355+1e0*-1e6
    df[2][:nelx*nely]=-dv/n/0.355

#
#   f[0] = v#/2121.
#   df[0][:nelx*nely] = dv#/2121.
#   df[0][-2:] = 0.
#
#   some plots
#
    itr=aux[-1]
    plot(nelx,nely,x_k,b_k,Gx,Gy,H,Hs,ovrh,itr)
    aux[-1]=aux[-1]+1
#
    return f, df
#
