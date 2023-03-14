#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'small',
          'figure.figsize': (25, 6),
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small',
         'font.family':'serif'}
pylab.rcParams.update(params)
#
def cost2d_simu(x_k, b_k, Gx, Gy, H, Hs):
#
#   filter
#
    x_f=np.zeros_like(x_k)
#   x_f[:]=np.asarray(H*x_k[np.newaxis].T)[:,0]
    x_f[:]=np.asarray(H*x_k[np.newaxis].T/Hs)[:,0]
#
#   get spatial gradients
#
    srf_gdx = np.asarray(Gx*x_f[np.newaxis].T)[:,0]
    srf_gdy = np.asarray(Gy*x_f[np.newaxis].T)[:,0]
    nab=np.array([srf_gdx,srf_gdy])
    eps=1e-12
    mag= np.sqrt(np.power(srf_gdx,2.) + np.power(srf_gdy,2.) + eps )  #np.linalg.norm(nab, axis=0)+eps
    bag=np.linalg.norm(b_k, axis=0)#+eps
#   mag=np.sqrt(np.power(srf_gdx,2.)+np.power(srf_gdy,2.) + eps)
#   bag=np.sqrt(np.power(b_k[0],2.)+np.power(b_k[1],2.) + eps)
#
#   gradient projeted along the build direction, normalised
#
    nng=1.4*np.dot(nab.T,b_k)/bag
#   nng=np.dot(nab.T,b_k)/mag/bag
    dnngdx=np.dot(nab.T,np.array([1.,0.]))/(mag)
    dnngdy=np.dot(nab.T,np.array([0.,1.]))/(mag)
    ang=np.rad2deg(np.arccos(nng))
#
#   calculate cost and derivative to build direction parameters
#
    inp = 1./3.
#
    [C, dC] = cost(nng)
    f = np.sum(C*x_f**inp)
#   f = np.sum(C*mag**2.*x_f**inp)
    dCdx = np.sum(1.4*dC*(srf_gdx/bag-b_k[0]*np.dot(nab.T,b_k)/np.power(bag,3.))*x_f**inp)
    dCdy = np.sum(1.4*dC*(srf_gdy/bag-b_k[1]*np.dot(nab.T,b_k)/np.power(bag,3.))*x_f**inp)
#   dCdx = np.sum(dC*(srf_gdx/bag-b_k[0]*np.dot(nab.T,b_k)/np.power(bag,3.))*x_f**inp)
#   dCdy = np.sum(dC*(srf_gdy/bag-b_k[1]*np.dot(nab.T,b_k)/np.power(bag,3.))*x_f**inp)
#   dCdx = np.sum(dC*(srf_gdx/bag-b_k[0]*np.dot(nab.T,b_k)/np.power(bag,3.))/mag*x_f**inp)
#   dCdy = np.sum(dC*(srf_gdy/bag-b_k[1]*np.dot(nab.T,b_k)/np.power(bag,3.))/mag*x_f**inp)
#
    dfdb = np.zeros(2,dtype=float)
    dfdb[0] = dCdx
    dfdb[1] = dCdy
#
#   total surface area and its derivative
#
    s = np.sum(mag)
    tmp = H/Hs # filter matrix (not transposed)
#
    rhs=srf_gdx[np.newaxis]/mag
    dmagdnp1 = np.asarray(rhs*Gx*tmp)[0,:]
    rhs=srf_gdy[np.newaxis]/mag
    dmagdnp2 = np.asarray(rhs*Gy*tmp)[0,:]
#   dmagdnp1 = np.asarray(tmp*Gx.T*(srf_gdx/mag)[np.newaxis].T)[:,0]
#   dmagdnp2 = np.asarray(tmp*Gy.T*(srf_gdy/mag)[np.newaxis].T)[:,0]
#
    dsdx = dmagdnp1 + dmagdnp2
#
#   surface cost and its derivative
#
#   [C,dC] = cost(nng)
#   t = np.sum(C*mag)
#
    rhs=dC*(1.4*b_k[0]/bag)*x_f**inp
#   rhs=dC*(b_k[0]/bag - np.dot(nab.T,b_k)*nab[0]/np.power(mag,2.)/bag)*mag*x_f**inp
#   rhs=dC*(b_k[0]/mag/bag - np.dot(nab.T,b_k)*nab[0]/np.power(mag,3.)/bag)*x_f**inp
    dcosdnx = np.asarray(rhs*Gx*tmp)[0,:]
    rhs=dC*(1.4*b_k[1]/bag)*x_f**inp
#   rhs=dC*(b_k[1]/bag - np.dot(nab.T,b_k)*nab[1]/np.power(mag,2.)/bag)*mag*x_f**inp
#   rhs=dC*(b_k[1]/mag/bag - np.dot(nab.T,b_k)*nab[1]/np.power(mag,3.)/bag)*x_f**inp
    dcosdny = np.asarray(rhs*Gy*tmp)[0,:]
#
    rhs=C*srf_gdx[np.newaxis]*2.*x_f**inp
    dmagdnp1 = np.asarray(rhs*Gx*tmp)[0,:]
    rhs=C*srf_gdy[np.newaxis]*2.*x_f**inp
    dmagdnp2 = np.asarray(rhs*Gy*tmp)[0,:]
#
    rhs=C*inp*x_f**(inp-1.)#*mag**2.
    dmagdxf = np.asarray(rhs*tmp)[0,:]
#
#   dfdx = dmagdnp1+dmagdnp2+dmagdxf
#   dfdx = dcosdnx+dcosdny+dmagdxf+dmagdnp1+dmagdnp2
    dfdx = dcosdnx+dcosdny+dmagdxf
#
#   measure of interior
#
    tmp = x_f*(1e-3/(mag+1e-3))
    i = np.sum(tmp)
    didx = np.zeros_like(dfdx)
#
    tmp = H/Hs # filter matrix (not transposed)
    rhs = 1e-3/(1e-3+mag) 
    rhs1= -1e-3*x_f/(1e-3 + mag)**2. *srf_gdx[np.newaxis]/mag 
    rhs2= -1e-3*x_f/(1e-3 + mag)**2. *srf_gdy[np.newaxis]/mag
    dintdx = np.asarray(rhs*tmp)[0,:]
    dint1dx = np.asarray(rhs1*Gx*tmp)[0,:]
    dint2dx = np.asarray(rhs2*Gy*tmp)[0,:]
#
    didx = dintdx + dint1dx + dint2dx
#
    return f, dfdb, dfdx, i, didx
#
def cost(nng):
#
#   domain is -1 to 1, with 1 being most severely overhanging
#
#   p=np.array([ 1.50000000e-01, -7.35931288e-04,  7.75000000e-01,  4.25735931e-01, -3.50000000e-01])
    p=np.array([ 0.15,        0.10611505,  0.9261101,   0.31888495, -0.5011101 ])
#
#   p=np.array([1., 1.])
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
    nng=np.dot(nab.T,b_k)/bag
#   nng=np.dot(nab.T,b_k)/mag/bag
    ang=np.rad2deg(np.arccos(nng))
#
#   Subplots
#
    fig,axs = plt.subplots(2,4)
#   major_xticks = np.array(major_xticks)
#   major_yticks = np.array(major_yticks)
    for j in range(2):
        for i in range(4):
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
    ovr=np.where(nng/mag > ovrh, 1., 0.)*mag
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
    [cst,_]=cost(nng)*x_f**(1./3.)
#
    axs[1,2].set_title(r'Cost (positive scalar) per surface area density e.g. $ C[\cos(\theta)] = \ldots $')
    cmap = plt.cm.rainbow
#   norm = colors.BoundaryNorm(np.arange(0, 180, 5), cmap.N)
    im = axs[1,2].imshow(cst.reshape((nely,nelx)), origin='lower',\
    extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap)#, vmin=0,vmax=1)
#   t = np.linspace(0., 1., num=8, endpoint=True)
    plt.colorbar(im, ax=axs[1,2])#, ticks=t)
#
#   interior
#
    rmin=2.
    tmp=np.power(x_f,1.)*(1e-3/(mag+1e-3))
    axs[0,3].set_title(r'Interior $ \mu = \rho {\epsilon}/ \left( \| \nabla \rho \| + \epsilon \right) $ %f (%f)'%(np.sum(tmp),(2*nelx/4-2*rmin)*(2*nely/8-2*rmin)))
    im = axs[0,3].imshow(tmp.reshape((nely,nelx)), cmap='gray_r', origin='lower',\
    interpolation='none',extent=(-nelx/2,nelx/2,-nely/2,nely/2),vmin=np.amin(tmp),vmax=np.amax(tmp))
    plt.colorbar(im, ax=axs[0,3])
#
#   Cost of surface
#
#    [cst,_]=cost(nng)*mag
#
    pen=2.
    icst= (2.-4.*pen)*tmp**2. + (4.*pen-1.)*tmp
    axs[1,3].set_title(r'Cost (positive scalar) per volume density e.g. $ C_V[\mu] = \ldots $')
    cmap = plt.cm.rainbow
#   norm = colors.BoundaryNorm(np.arange(0, 180, 5), cmap.N)
    im = axs[1,3].imshow(icst.reshape((nely,nelx)), origin='lower',\
    extent=(-nelx/2,nelx/2,-nely/2,nely/2), cmap=cmap)#, vmin=0,vmax=1)
#   t = np.linspace(0., 1., num=8, endpoint=True)
    plt.colorbar(im, ax=axs[1,3])#, ticks=t)
#
#   stop
#
    plt.figtext(0.015, 0.025, r'Sum E-norm of spatial gradient (surface area) $ \int_V \| \nabla \rho \| d V = $ %.1f'%(np.sum(mag)), fontsize=13)
    plt.figtext(0.27, 0.025, r'Sum E-norm overhanging spatial gradient (overhanging s.area) $ \int_O \| \nabla \rho \| d V = $  %.1f'%(np.sum(ovr)), fontsize=13)
    plt.figtext(0.6, 0.025, r'Cost, function of angle w.r.t. b-direct. ... $ \int_V C[\theta] \| \nabla \rho \| d V =  $  %.1f'%(np.sum(cst)), fontsize=13)
#
    plt.figtext(0.85, 0.025, r'Volume cost ... $ \int_V C_V[\mu] d V =  $  %.1f'%(np.sum(icst)), fontsize=13)
#
    fig.tight_layout(pad=4)
    plt.savefig('./png/topo_%03d.png'%k, dpi=99)#, dpi=199)
    plt.close('all')
#
