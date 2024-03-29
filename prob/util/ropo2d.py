#
import os
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
import cvxopt; import cvxopt.cholmod
#
def topo2d_init(nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,fixed,force,pen,qen,muc,Emin,Emax,gv,g):
#
    # dofs
    nelm=nelx*nely
    ndof = 2*(nelx+1)*(nely+1)

    # FE: Build the index vectors for the for coo matrix format.
    KE=lk()
    edofMat=np.zeros((nelx*nely,8),dtype=int)
    for ely in range(nely):
        for elx in range(nelx):
            el = ely*nelx+elx
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
    for j in range(fely):
        for i in range(felx):
            row=j*felx+i
            kk1=int(np.maximum(j-(np.ceil(rmin)-1),0))
            kk2=int(np.minimum(j+np.ceil(rmin),fely))
            ll1=int(np.maximum(i-(np.ceil(rmin)-1),0))
            ll2=int(np.minimum(i+np.ceil(rmin),felx))
            for k in range(kk1,kk2):
                for l in range(ll1,ll2):
                    col=k*felx+l
                    fac=rmin-np.sqrt(((i-l)*(i-l)+(j-k)*(j-k)))
                    iH[cc]=row
                    jH[cc]=col
                    sH[cc]=np.maximum(0.0,fac)
                    cc=cc+1
    # Finalize assembly and convert to csc format
    H=coo_matrix((sH,(iH,jH)),shape=(felx*fely,felx*fely)).tocsc()    
    Hs=H.sum(1)
#
    l=2.
    ngrad=int((nely-1)*(nelx-1)*2)
#
    iG = np.zeros(ngrad)
#
    jGx = np.zeros(ngrad)
    sGx = np.zeros(ngrad)
    jGy = np.zeros(ngrad)
    sGy = np.zeros(ngrad)
#
    cc=0
    for j in range(1,nely-1):
        for i in range(1,nelx-1):
#
            c = j*nelx+i
            cr = j*nelx+i+1
            cl = j*nelx+i-1
            cu = (j+1)*nelx+i
            cd = (j-1)*nelx+i
#
            iG[cc]=c
            iG[cc+1]=c
#
            jGx[cc]=cr
            jGx[cc+1]=cl
            jGy[cc]=cu
            jGy[cc+1]=cd
#
            sGx[cc]=1./l
            sGx[cc+1]=-1./l
#
            sGy[cc]=1./l
            sGy[cc+1]=-1./l
#
            cc=cc+2
#
    Gx=coo_matrix((sGx,(iG,jGx)),shape=(nelx*nely,nelx*nely)).tocsc()
    Gy=coo_matrix((sGy,(iG,jGy)),shape=(nelx*nely,nelx*nely)).tocsc()
#
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
    aux=[nelx,nely,v_l,v_0,v_u,ft,rmin,felx,fely,xPadd,pen,qen,muc,Emin,Emax,gv,\
        ndof,KE,H,Hs,Gx,Gy,iK,jK,edofMat,fixed,free,f,u,0]
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
        ndof,KE,H,Hs,Gx,Gy,iK,jK,edofMat,fixed,free,f,u,itr]=aux
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
    # Setup and solve FE problem
    sK=((KE.flatten()[np.newaxis]).T*(Emin+xPena(muc,pen,xPhys)*(Emax-Emin))).flatten(order='F')
    K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K,fixed,fixed).tocoo()
    # Set self-weight load
    f_tmp=np.zeros(ndof)
    np.add.at(f_tmp, edofMat[:, 1::2].flatten(), np.kron(xPena(0.,qen,xPhys), gv * np.ones(4)/4. ))
    f_apl = f.copy()
    f_apl[:,0] += f_tmp
    # Solve system 
    K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
    B = cvxopt.matrix(f_apl[free,0])
    cvxopt.cholmod.linsolve(K,B)
    u[free,0]=np.array(B)[:,0]

    # Objective and sensitivity
    ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE)*u[edofMat].reshape(nelx*nely,8) ).sum(1)
    obj = ( (Emin+xPena(muc,pen,xPhys)*(Emax-Emin))*ce ).sum()
    tmp = -dxPena(muc,pen,xPhys)*(Emax-Emin)*ce
    tmp += 2. * gv * u[edofMat[:,1],0] * 1./4. * dxPena(0.,qen,xPhys)#qen*xPhys**(qen-1.)
    tmp += 2. * gv * u[edofMat[:,3],0] * 1./4. * dxPena(0.,qen,xPhys)#qen*xPhys**(qen-1.)
    tmp += 2. * gv * u[edofMat[:,5],0] * 1./4. * dxPena(0.,qen,xPhys)#qen*xPhys**(qen-1.)
    tmp += 2. * gv * u[edofMat[:,7],0] * 1./4. * dxPena(0.,qen,xPhys)#qen*xPhys**(qen-1.)
    dc[:] = tmp
    dv[:] = np.ones(nely*nelx)

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
#element stiffness matrix
def lk():
    E=1
    nu=0.3
    k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
    KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
    [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
    [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
    [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
    [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
    [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
    [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
    [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
    return (KE)
#
def deleterowcol(A, delrow, delcol):
    # Assumes that matrix is in symmetric csc form !
    m = A.shape[0]
    keep = np.delete (np.arange(0, m), delrow)
    A = A[keep, :]
    keep = np.delete (np.arange(0, m), delcol)
    A = A[:, keep]
    return A    
