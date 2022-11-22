#
import numpy as np
import importlib.util
from subs.t2cplx import t2c as subs
#from subs.t2duel import t2d as subs
#from prob.util.orien import orien_init
#from prob.util.orien import orien_simu
#from prob.util.orien import orien_outp
#
spec=importlib.util.spec_from_file_location("orien", "/home/dirk/RECIPE/knap/orien.py")
orien = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orien)
#
def init(g):
#
#   Run init from RECIPE module
    [ply,nrm,aea,fln]=orien.orien_init("/home/dirk/RECIPE/knap/stl/dunny.stl")
#
    t=len(aea)
#
#   Number of variables is 4 quaternion coefficients (put them at the end), 
#   and then an additional 2 per triangle (one for down indicator, and one for overhanging indicator)
#   the length of the array of triangle areas is thus used
    n = 2*t + 4
#
#   The number of constraints is 2 per triangle; one to set the down indicator, and one to set the 
#   overhanging indicator; plus one more to constrain the quaternion coefficients to length 1
    m = 2*t + 1
#
#   Indicator variables, followed by 4 quaternion coefficients
    x_l = np.ones(n,dtype=float)*0.0
    x_u = np.ones(n,dtype=float)
    x_k = np.ones(n,dtype=float)*0.5
    x_l[-4:] = -1.
#
    x_k[-4] = 0.
    x_k[-3] = 0.
    x_k[-2] = 0.
    x_k[-1] = 0.
#
    tmp=np.load('glob_48.npz')
    x_k[:] = tmp['x_i']
#
#   Last constraint is the quaternion norm (equality)
    c_s = np.ones(m,dtype=int)
    c_s[-1] = 0
#
#   Pack data returned from init into aux
    aux=[ply,nrm,aea,t,fln,0]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def simu(n,m,x,aux,g):
#
    f = np.zeros((m + 1), dtype=float)
    df = [np.zeros(n,dtype=float)]
#
    q = x[-4:].copy()
    [ply,nrm,aea,t,fln,k]=aux
    [rrm,dndr,dndi,dndj,dndk]=orien.orien_simu(q,nrm)
#
#   x = [y,z,q]
#
    f[0] = np.sum(x[:t]*x[t:2*t]*aea)/np.sum(aea)
#
    tmp = rrm.copy()
    tmp[:,0]=-2.*tmp[:,0]; tmp[:,1]=-2.*tmp[:,1]; tmp[:,2]= 2.*tmp[:,2]
    dfdr = tmp*dndr
    dfdi = tmp*dndi
    dfdj = tmp*dndj
    dfdk = tmp*dndk
    p = 9.
    for j in range(t):
        df[0][j] = x[t+j]*aea[j]/np.sum(aea)
        df[0][t+j] = x[j]*aea[j]/np.sum(aea)
        f[j+1] = -rrm[j,2] - x[j]**p
        df.append((j,j,-p*x[j]**(p-1.)))
        df.append((j,n-4,-dndr[j,2]))
        df.append((j,n-3,-dndi[j,2]))
        df.append((j,n-2,-dndj[j,2]))
        df.append((j,n-1,-dndk[j,2]))
        f[t+j+1] = rrm[j,2]**2. - rrm[j,1]**2. - rrm[j,0]**2. - x[t+j]**p
        df.append((t+j,t+j,-p*x[t+j]**(p-1.)))
        df.append((t+j,n-4, 2.*rrm[j,2]*dndr[j,2]-2.*rrm[j,1]*dndr[j,1]-2.*rrm[j,0]*dndr[j,0]))
        df.append((t+j,n-3, 2.*rrm[j,2]*dndi[j,2]-2.*rrm[j,1]*dndi[j,1]-2.*rrm[j,0]*dndi[j,0]))
        df.append((t+j,n-2, 2.*rrm[j,2]*dndj[j,2]-2.*rrm[j,1]*dndj[j,1]-2.*rrm[j,0]*dndj[j,0]))
        df.append((t+j,n-1, 2.*rrm[j,2]*dndk[j,2]-2.*rrm[j,1]*dndk[j,1]-2.*rrm[j,0]*dndk[j,0]))
#
#   unit quaternion constraint and derivative terms
    f[-1] = q[0]**2. + q[1]**2. + q[2]**2. + q[3]**2. - 1.0
    df.append((m-1,n-4,2.*q[0]))
    df.append((m-1,n-3,2.*q[1]))
    df.append((m-1,n-2,2.*q[2]))
    df.append((m-1,n-1,2.*q[3]))
#
    kt=aux[-1]
    kt=kt+1
    aux[-1]=kt
    orien.orien_outp(fln,ply,x,t,kt)
#
    over=np.where(rrm[:,2]<0,1.0,0)*np.where(-rrm[:,2]>np.linalg.norm(rrm[:,:2],axis=1),1.0,0.)
    fd=np.sum(over*aea)/np.sum(aea)
    print('F: %14.7e Q: (%4.1f,%4.1f,%4.1f,%4.1f)'%(fd,q[0],q[1],q[2],q[3]))
#
    return f, df
#
def apar(n):
#
    mov=1.0*np.ones(n)
    asf=[0.7,1.1]
#
    enf='None'
#       
    kmx=1000
    cnv=[1e-3,1e-3,1e4,1e-6,1e-3]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[np.ones_like(x_k)]
    if k > 0:
        sph = f_1 - f_k
        sph[0]=sph[0]-np.dot(df_k[0],x_1-x_k)
        for df in df_k[1:]: sph[df[0]+1]=sph[df[0]+1]-df[2]*(x_1[df[1]]-x_k[df[1]])
        sph=2.*sph/np.maximum(np.linalg.norm(x_1-x_k)**2.,1e-6)
        c_x[0]=c_x[0]*sph[0]
        for j in range(len(f_k)-1):
            for i in range(len(x_k)):
                c_x.append((j,i,sph[j+1]))
    else:
        for j in range(len(f_k)-1):
            for i in range(len(x_k)):
                c_x.append((j,i,1.))
#
    d_l= np.maximum(x_l, x_k-mov*(x_u-x_l))
    d_u= np.minimum(x_u, x_k+mov*(x_u-x_l))
#
    L=L_k
    U=U_k
#
    return c_x,mov,L,U,d_l,d_u
#
