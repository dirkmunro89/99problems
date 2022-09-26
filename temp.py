#
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt; import cvxopt.cholmod
#
# specify subsolver here
#
from subs.condual import con as subs
#
# specify problem and algorithmic parameters here
#
def apar():
#   
    mov=0.5
    asf=[0.7,1.1]
#       
    kmx=200
    cnv=[1e-6,1e-6]
#       
    return mov, asf, kmx, cnv
#
def caml(k, x_k, dg, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(dg)
#
    if k<=1:
        L = x_k-mov*(x_u-x_l)
        U = x_k+mov*(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - asf[0]*(x_1 - L_k), x_k - asf[1]*(x_1 - L_k))
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + asf[0]*(U_k - x_1), x_k + asf[1]*(U_k - x_1))
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,L,U,d_l,d_u
#
def init():
#
    nelx=20
    nely=20
    volfrac_lo=0.01
    volfrac_0=0.8
    volfrac_up=0.8
    rmin=1.1
    penal=2.0
    ft=1 # ft==0 -> sens, ft==1 -> dens
 
    n = nelx*nely; m = 2
    x_l = np.ones(n,dtype=float)*0.2
    x_u = np.ones(n,dtype=float)
    x_k = volfrac_0*np.ones(n,dtype=float)
 
    # Max and min stiffness
    Emin=0e-9; Emax=1.0
    gv=-9.81/1e3

    # dofs
    ndof = 2*(nelx+1)*(nely+1)

    # FE: Build the index vectors for the for coo matrix format.
    KE=lk()
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
    nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
    iH = np.zeros(nfilter)
    jH = np.zeros(nfilter)
    sH = np.zeros(nfilter)
    cc=0
    for i in range(nelx):
        for j in range(nely):
            row=i*nely+j
            kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
            kk2=int(np.minimum(i+np.ceil(rmin),nelx))
            ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
            ll2=int(np.minimum(j+np.ceil(rmin),nely))
            for k in range(kk1,kk2):
                for l in range(ll1,ll2):
                    col=k*nely+l
                    fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
                    iH[cc]=row
                    jH[cc]=col
                    sH[cc]=np.maximum(0.0,fac)
                    cc=cc+1
    # Finalize assembly and convert to csc format
    H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()    
    Hs=H.sum(1)

    # BC's and support
    dofs=np.arange(2*(nelx+1)*(nely+1))
    ndofy=2*(nely+1)
    fixed = np.union1d(dofs[0:ndofy:2],np.array([ndof - 2 , ndof - 1]))
    free=np.setdiff1d(dofs,fixed)

    # Solution and RHS vectors
    f=np.zeros((ndof,1))
    u=np.zeros((ndof,1))

    # Set load
    f[1,0]=0

    xPhys=np.ones(n,dtype=float)*volfrac_0

    # Initialize plot and plot the initial design
    plt.ion()
    fig,ax = plt.subplots()
    im = ax.imshow(-xPhys.reshape((nelx,nely)).T, cmap='gray',\
    interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
    fig.show()
#
    aux=[nelx,nely,volfrac_lo,volfrac_0,volfrac_up,rmin,penal,ft,Emax,Emin,gv,\
        ndof,KE,H,Hs,iK,jK,edofMat,fixed,free,f,u,im,fig]
#
    return n,m,x_l,x_u,x_k,aux
#
def simu(n,m,x,aux):
#
    g = np.zeros((m + 1), dtype=float)
    dg = np.zeros((m + 1, n), dtype=float)
#
    ce=np.zeros(n,dtype=float)
    dc=np.zeros(n,dtype=float)
    xPhys=np.zeros(n,dtype=float)
    dv=np.ones(n,dtype=float)
#
    [nelx,nely,volfrac_lo,volfrac_0,volfrac_up,rmin,penal,ft,Emax,Emin,gv,\
        ndof,KE,H,Hs,iK,jK,edofMat,fixed,free,f,u,im,fig]=aux
#
    # Filter design variables
    if ft==0:   xPhys[:]=x
    elif ft==1:    xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]
#
    # Plot to screen and save
    im.set_array(-xPhys.reshape((nelx,nely)).T)
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.savefig('topo.eps')
#
#   muc=0.25
#   xPena = np.where(xPhys < muc, xPhys*muc**(penal-1e0), xPhys**penal   )
#   dxPena = np.where(xPhys < muc, muc**(penal-1e0), penal*xPhys**(penal-1)   )
#
    # Setup and solve FE problem
#   sK=((KE.flatten()[np.newaxis]).T*((xPena)*(Emax-Emin))).flatten(order='F')
#   sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPena)*(Emax-Emin))).flatten(order='F')
#   sK=((KE.flatten()[np.newaxis]).T*((xPena)*(Emax))).flatten(order='F')
    sK=((KE.flatten()[np.newaxis]).T*((xPhys)**penal*(Emax))).flatten(order='F')
    K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
    # Remove constrained dofs from matrix and convert to coo
    K = deleterowcol(K,fixed,fixed).tocoo()
    # Set self-weight load
    f_tmp=np.zeros(ndof)
    np.add.at(f_tmp, edofMat[:, 1::2].flatten(), np.kron(xPhys, gv * np.ones(4)/4. ))
    f_apl = f.copy()
    f_apl[:,0] += f_tmp
    # Solve system 
    K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
    B = cvxopt.matrix(f_apl[free,0])
    cvxopt.cholmod.linsolve(K,B)
    u[free,0]=np.array(B)[:,0]

    # Objective and sensitivity
    ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE)*u[edofMat].reshape(nelx*nely,8) ).sum(1)
#   obj = ( (xPena*(Emax))*ce ).sum()
#   obj = ( (Emin+xPena*(Emax-Emin))*ce ).sum()
    obj = ( (xPhys**penal*(Emax))*ce ).sum()
#   dc[:] = (-dxPena*(Emax))*ce
#   dc[:] = (-dxPena*(Emax-Emin))*ce
    dc[:] = (-penal*xPhys**(penal-1)*(Emax))*ce
    dc[:] += 2. * gv * u[edofMat[:,1],0] * 1./4.; dc[:] += 2. * gv * u[edofMat[:,3],0] * 1./4.
    dc[:] += 2. * gv * u[edofMat[:,5],0] * 1./4.; dc[:] += 2. * gv * u[edofMat[:,7],0] * 1./4.
    dv[:] = np.ones(nely*nelx)
    # Sensitivity filtering:
    if ft==0:
        dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
    elif ft==1:
        dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
        dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]
#
    g[0]=obj
    g[1]=np.sum(x)/n-volfrac_up
    g[2]=-np.sum(x)/n+volfrac_lo
#
#   print(np.sum(x)/n)
#
    dg[0][:] = dc
    dg[1][:] = dv/n
    dg[2][:] = -dv/n
#
    return g, dg
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
#
