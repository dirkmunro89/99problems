#
import os
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt; import cvxopt.cholmod
#
from prob.util.para2d import para2d
#
def topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fixed,force,pen,qen,muc,Emin,Emax,gv,g):
#
    # dofs
    nelm=nelx*nely
    ndof = 2*(nelx+1)*(nely+1)

    # FE: Build the index vectors for the for coo matrix format.
    E=1
    nu=0.3
    l=1.
    t=1.
    DE=ld(l)
    SE=ls(E,nu,l,t)
    NE=ln(l,t)
#
    edofMat=np.zeros((nelx*nely,8),dtype=int)
    for elx in range(nelx):
        for ely in range(nely):
            el = ely+elx*nely
            n1=(nely+1)*elx+ely
            n2=(nely+1)*(elx+1)+ely
            edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
    # Construct the index pointers for the coo format
    iK = np.kron(edofMat,np.ones((8,1))).flatten()
    jK = np.kron(edofMat,np.ones((1,8))).flatten()    

    # Filter: Build (and assemble) the index+data vectors for the coo matrix format
    nfilter=int(felx*fely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc=0
    for i in range(felx):
        for j in range(fely):
            row=i*fely+j
            kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
            kk2=int(np.minimum(i+np.ceil(rmin),felx))
            ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
            ll2=int(np.minimum(j+np.ceil(rmin),fely))
            for k in range(kk1,kk2):
                for l in range(ll1,ll2):
                    col=k*fely+l
                    fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
                    iH[cc]=row
                    jH[cc]=col
                    sH[cc]=np.maximum(0.0,fac)
                    cc=cc+1
    # Finalize assembly and convert to csc format
    H=coo_matrix((sH,(iH,jH)),shape=(felx*fely,felx*fely)).tocsc()    
    Hs=H.sum(1)

    # BC's and support
    dofs=np.arange(ndof)
#   fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
    free=np.setdiff1d(dofs,fixed)

    # Solution and RHS vectors
    f=np.zeros((ndof,1))
    u=np.zeros((ndof,1))

    # Set load
    f[[i[0] for i in force],0]=[i[1] for i in force]

#
    if g <= 0 and g > -2:
        plt.ion(); fig,ax = plt.subplots()
        im = ax.imshow(-np.zeros(nelm).reshape((nelx,nely)).T, cmap='gray',\
        interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
    else:
        fig=0; im=0
#
    aux=[nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,pen,qen,muc,Emin,Emax,gv,\
        ndof,SE,DE,NE,H,Hs,iK,jK,edofMat,fixed,free,f,u,fig,im]
#
    return aux
#
def xPena(muc,penal,xPhys):
    return muc*xPhys + (1-muc)*xPhys**penal
def dxPena(muc,penal,xPhys):
    return muc + penal*(1-muc)*xPhys**(penal-1.)
#
def topo2d_simu(n,m,x,aux,vis):
#
    ce=np.zeros(n,dtype=float)
    dc=np.zeros(n,dtype=float)
    xPhys=np.zeros(n,dtype=float)
    dv=np.ones(n,dtype=float)
#
    [nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,pen,qen,muc,Emin,Emax,gv,\
        ndof,SE,DE,NE,H,Hs,iK,jK,edofMat,fixed,free,f,u,fig,im]=aux
#
    if np.count_nonzero(xPadd <= -1) != len(x):
        print('Filter padding error')
        stop
#
    pad=np.argwhere(xPadd <= -1).flatten()
#
# Filter design variables
    if ft==0:   xPhys[:]=x
    elif ft==1: 
        tmp=xPadd.copy(); tmp[pad]=x
        xPhys[:]=np.asarray(H*tmp[np.newaxis].T/Hs)[:,0][pad] ###
#
    if fig and im:
        im.set_array(-x.reshape((nelx,nely)).T)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.savefig('topology.eps')
        tmp=np.flip(x.reshape((nelx,nely)).T,0)
        tmp2=np.append(np.flip(tmp,1),tmp,axis=1)
        np.savetxt("topology.dat",tmp2,fmt='%14.7f')
    else:
        if vis > 0:
            fig,ax = plt.subplots()
            im = ax.imshow(-x.reshape((nelx,nely)).T, cmap='gray',\
            interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
            plt.savefig('topo_%d.eps'%vis)
            plt.close()
#
    nelm=nelx*nely
    sK=np.zeros((nelm*64),dtype=float)
    fal = f.copy()
    bf=np.array([[0.],[gv],[0.],[gv],[0.],[gv],[0.],[gv]])
#
#   'assembly', slow if done in this way, of course, but wanting to maintain a vanilla 
#   loop-based implenentation for going to C
#
    for e in range(nelm):
        SP=SE*(Emin+xPena(muc,pen,xPhys[e])*(Emax-Emin))
        KE=np.dot(DE.T,np.dot(SP,DE))
        sK[64*e:64*(e+1)]=KE.flatten()
        fal[edofMat[e]] = fal[edofMat[e]] + np.dot(NE,bf)*xPena(0.,qen,xPhys[e])
#
    K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K,fixed,fixed).tocoo()
    # Solve system 
    K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
    B = cvxopt.matrix(fal[free,0])
    cvxopt.cholmod.linsolve(K,B)
    u[free,0]=np.array(B)[:,0]
#
#   'post processing' (obj and sens), also slow if done in this way, of course; see comment above
#
    obj=0.
    tmp=np.zeros(n,dtype=float)
    for e in range(nelm):
        SP=SE*(Emin+xPena(muc,pen,xPhys[e])*(Emax-Emin))
        eps=np.dot(DE,u[edofMat[e]])
        ele=np.dot(eps.T,np.dot(SP,eps))[0][0]
        ele0=np.dot(eps.T,np.dot(SE,eps))[0][0] # can also be done without redoing this
        obj=obj+ele
        tmp[e]=-ele0*dxPena(muc,pen,xPhys[e])
        tmp[e]=tmp[e] + 2.*np.dot(u[edofMat[e]].T,np.dot(NE,bf))
    dc[:] = tmp
    dv[:] = np.ones(nely*nelx)

#   c=0
#   while os.path.isfile('u_vec_%d.dat'%c): c=c+1
#   np.savetxt('u_vec_%d.dat'%c, u[edofMat[:,[1,3,5,7]],0] )
#   para2d(nelx,nely,x,u,gv,tmp)

    # Sensitivity filtering:
    if ft==0:
        dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
    elif ft==1:
        tmp=np.zeros_like(xPadd)
        tmp[pad]=dc; dc[:] = np.asarray(H*(tmp[np.newaxis].T/Hs))[:,0][pad]
        tmp[pad]=dv; dv[:] = np.asarray(H*(tmp[np.newaxis].T/Hs))[:,0][pad]
#
    v=np.sum(xPhys)
#
    return obj,dc,v,dv
#
def deleterowcol(A, delrow, delcol):
    # Assumes that matrix is in symmetric csc form !
    m = A.shape[0]
    keep = np.delete (np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete (np.arange(0, m), delcol)
    A = A[:, keep]
    return A    
#
def ls(E,nu,l,t):
#
#
    s=E*l**2.*t/4./(1.-nu**2.)*np.array([1., nu, (1.-nu)/2.])
    SE = np.zeros((12,12),dtype=float)
#
    for i in range(4):
        SE[i*3+0][i*3+0]=s[0]; SE[i*3+0][i*3+1]=s[1]
        SE[i*3+1][i*3+0]=s[1]; SE[i*3+1][i*3+1]=s[0]
        SE[i*3+2][i*3+2]=s[2]
#
    return SE
#
def ld(l):
    d=np.array([-np.sqrt(3e0)-1e0,-np.sqrt(3e0)+1e0,np.sqrt(3e0)-1e0,np.sqrt(3e0)+1e0])/np.sqrt(3)/2./l
    DE = np.array([
    [d[0],   0., d[3],   0., d[2],   0., d[1],   0.],
    [  0., d[0],   0., d[1],   0., d[2],   0., d[3]],
    [d[0], d[0], d[1], d[3], d[2], d[2], d[3], d[1]],
#
    [d[0],   0., d[3],   0., d[2],   0., d[1],   0.],
    [  0., d[1],   0., d[0],   0., d[3],   0., d[2]],
    [d[1], d[0], d[0], d[3], d[3], d[2], d[2], d[1]],
#
    [d[1],   0., d[2],   0., d[3],   0., d[0],   0.],
    [  0., d[1],   0., d[0],   0., d[3],   0., d[2]],
    [d[1], d[1], d[0], d[2], d[3], d[3], d[2], d[0]],
#
    [d[1],   0., d[2],   0., d[3],   0., d[0],   0.],
    [  0., d[0],   0., d[1],   0., d[2],   0., d[3]],
    [d[0], d[1], d[1], d[2], d[2], d[3], d[3], d[0]],
    ]);
#
    return (DE)
#
def ln(l,t):
    n=np.array([(np.sqrt(3)+1)**2.,(np.sqrt(3)-1.)*(np.sqrt(3)+1.), (np.sqrt(3)-1.)**2.])*l**2.*t/48.
    NE = np.array([
    [n[0],   0., n[1],   0., n[2],   0., n[1],   0.],
    [  0., n[0],   0., n[1],   0., n[2],   0., n[1]],
#
    [n[1],   0., n[0],   0., n[1],   0., n[2],   0.],
    [  0., n[1],   0., n[0],   0., n[1],   0., n[2]],
#
    [n[2],   0., n[1],   0., n[0],   0., n[1],   0.],
    [  0., n[2],   0., n[1],   0., n[0],   0., n[1]],
#
    [n[1],   0., n[2],   0., n[1],   0., n[0],   0.],
    [  0., n[1],   0., n[2],   0., n[1],   0., n[0]],
    ]).T;
#
    return (NE)
