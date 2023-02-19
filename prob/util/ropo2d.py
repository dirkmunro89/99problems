#
import os
import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
          'figure.figsize': (20, 11),
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small',
         'font.family':'serif'}
pylab.rcParams.update(params)
import cvxopt; import cvxopt.cholmod
#
#from prob.util.para2d import para2d
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
#   if g <= 0 and g > -2:
#       plt.ion(); fig,ax = plt.subplots()
#       im = ax.imshow(-np.zeros(nelm).reshape((nelx,nely)).T, cmap='gray_r',\
#       interpolation='none',norm=colors.Normalize(vmin=0,vmax=1), origin='lower')
#   else:
#       fig=0; im=0
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
#   if fig and im:
#       im.set_array(x.reshape((nely,nelx)))
#       fig.canvas.draw()
#       fig.canvas.flush_events()
#       plt.savefig('topology.eps')
#       tmp=np.flip(x.reshape((nely,nelx)),0)
#       tmp2=np.append(np.flip(tmp,1),tmp,axis=1)
#       np.savetxt("topology.dat",tmp2,fmt='%14.7f')
#   else:
#       if vis > 0:
#           fig,ax = plt.subplots()
#           im = ax.imshow(-x.reshape((nelx,nely)).T, cmap='gray',\
#           interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
#           plt.savefig('topo_%d.eps'%vis)
#           plt.close()
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
#   call cost stuff
#

#
#
#
    x_k=np.zeros(nelx*nely,dtype=float)
    x_k[:]=x[:nelx*nely]
    itr=aux[-1]
    ovrh = np.cos(np.pi/180.*45.)
    b_k=np.array([1.,0.])
    plot(nelx,nely,x_k,b_k,Gx,Gy,H,Hs,ovrh,itr)
    aux[-1]=aux[-1]+1
#
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
#
def plot(nelx,nely,x_k,b_k,Gx,Gy,H,Hs,ovrh,k):
#
    x_f=np.zeros_like(x_k)
    x_f[:]=np.asarray(H*x_k[np.newaxis].T/Hs)[:,0]
#
    srf_gdx = np.asarray(Gx*x_f[np.newaxis].T)[:,0]
    srf_gdy = np.asarray(Gy*x_f[np.newaxis].T)[:,0]
    nab=np.array([srf_gdx,srf_gdy])
    eps=1e-12
    mag=np.linalg.norm(nab, axis=0)+eps
    bag=np.linalg.norm(b_k, axis=0)+eps
    nng=np.dot(nab.T,b_k)/mag/bag
    ang=np.rad2deg(np.arccos(nng))
#
#   Subplots
#
    fig,axs = plt.subplots(2,3)
#   major_xticks = np.array(major_xticks)
#   major_yticks = np.array(major_yticks)
    for j in range(2):
        for i in range(3):
            axs[j,i].arrow(0,0,b_k[0]/bag*nelx/10.,b_k[1]/bag*nely/10,color='magenta',head_width=1.*nelx/20.)
            axs[j,i].grid()
#           if nelx < 11 and nely < 11:
#               axs[j,i].set_xticks(major_xticks)
#               axs[j,i].set_yticks(major_yticks)
#               axs[j,i].grid(which='both')
#               axs[j,i].grid(which='major', alpha=0.5)
#
#   Topology plot
#
    axs[0,0].set_title(r'Decision density field $\rho$ (possibly filtered)')
    im = axs[0,0].imshow(x_f.reshape((nely,nelx)), cmap='gray_r', origin='lower',\
    interpolation='none',extent=(-nelx/2,nelx/2,-nely/2,nely/2))#,vmin=np.amin(x_f),vmax=np.amax(x_f))
    plt.colorbar(im, ax=axs[0,0])
#
    axs[1,0].set_title(r'Norm of spatial gradient (surface area density) $\| \nabla \rho \|$ (Euclidean)')
    im = axs[1,0].imshow(mag.reshape((nely,nelx)), cmap='gray_r',origin='lower',\
    interpolation='none',extent=(-nelx/2,nelx/2,-nely/2,nely/2))
    plt.colorbar(im, ax=axs[1,0])
#
#   Ang. of variation plot
#
    axs[0,1].set_title(r'Angle with build direction $\theta = \cos^{-1}\left( b \cdot \nabla \rho /  \| \nabla \rho\|  \right)$ with $0 \leq \theta \leq 180^{\circ}$')
    cmap = plt.cm.rainbow
#   norm = colors.BoundaryNorm(np.arange(0, 180, 5), cmap.N)
    im = axs[0,1].imshow(ang.reshape((nely,nelx)), origin='lower',\
    extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap, vmin=0.,vmax=180.)
    t = np.linspace(0., 180, num=8, endpoint=True)
    plt.colorbar(im, ax=axs[0,1], ticks=t)
#
#   Overhang plot
#
    cmap = plt.cm.rainbow
    ovr=np.where(nng > ovrh, 1., 0.)*mag
#
    axs[1,1].set_title(r'Overhanging surface area density $\| \nabla \rho \|$ if $\cos(\theta) > \cos(45^{\circ} )$')
    im = axs[1,1].imshow(ovr.reshape((nely,nelx)), origin='lower',\
        extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap, vmin=0,vmax=1)
    plt.colorbar(im, ax=axs[1,1])
#
#   Ang. of variation plot
#
    axs[0,2].set_title(r'Normalised dot product with build direction $\cos(\theta) = { b \cdot \nabla \rho } / {  \| \nabla \rho \|  }  $')
    cmap = plt.cm.rainbow
#   norm = colors.BoundaryNorm(np.arange(0, 180, 5), cmap.N)
    im = axs[0,2].imshow(nng.reshape((nely,nelx)), origin='lower',\
    extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap)#, vmin=-1,vmax=1)
#   t = np.linspace(-1., 1., num=8, endpoint=True)
    plt.colorbar(im, ax=axs[0,2])#, ticks=t)
#
#   Cost of surface
#
    [cst,_]=cost(nng)*mag
#
    axs[1,2].set_title(r'Cost (positive scalar) per surface area density e.g. $ C[\cos(\theta)] = \ldots $')
    cmap = plt.cm.rainbow
#   norm = colors.BoundaryNorm(np.arange(0, 180, 5), cmap.N)
    im = axs[1,2].imshow(cst.reshape((nely,nelx)), origin='lower',\
    extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap)#, vmin=0,vmax=1)
#   t = np.linspace(0., 1., num=8, endpoint=True)
    plt.colorbar(im, ax=axs[1,2])#, ticks=t)
#
#   stop
#
    plt.figtext(0.015, 0.025, r'Sum E-norm of spatial gradient (surface area) $ \int_V \| \nabla \rho \| d V = $ %.1f'%(np.sum(mag)), fontsize=13)
    plt.figtext(0.315, 0.025, r'Sum E-norm overhanging spatial gradient (overhanging s.area) $ \int_O \| \nabla \rho \| d V = $  %.1f'%(np.sum(ovr)), fontsize=13)
    plt.figtext(0.7, 0.025, r'Cost, function of angle w.r.t. b-direct. ... $ \int_V C[\theta] \| \nabla \rho \| d V =  $  %.1f'%(np.sum(cst)), fontsize=13)
#
    fig.tight_layout(pad=4)
    plt.savefig('./png/topo_%03d.png'%k, dpi=99)#, dpi=199)
    plt.close('all')
#
def cost(nng):
#
#   domain is -1 to 1, with 1 being most severely overhanging
#
    p=np.array([ 1.50000000e-01, -7.35931288e-04,  7.75000000e-01,  4.25735931e-01, -3.50000000e-01])
#
    i=0; C=0.
    dC = 0.
    for c in p:
        C = C + c*np.power(nng,i)
        if i > 0:
            dC = dC + i*c*np.power(nng,i-1)
        i=i+1
#
    return C, dC
#
