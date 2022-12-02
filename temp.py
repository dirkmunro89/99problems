#
import numpy as np
import importlib.util
from subs.t2milx import t2c as subs
#from subs.t2duel import t2d as subs
#from prob.util.orien import orien_init
#from prob.util.orien import orien_simu
#from prob.util.orien import orien_outp
#
spec=importlib.util.spec_from_file_location("orien", "/home/dirk/RECIPE/knap/orien_suppt.py")
orien = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orien)
#
def init(g):
#
#   Run init from RECIPE module
    [ply,nrm,aea,cen,nog,fln]=orien.orien_init("/home/dirk/RECIPE/knap/stl/dunny.stl")
#
#   what to do here...? Get furthest point from COG, think of a sphere that goes around that point
#   to determine base-plate distance. maybe move cog so that that is zero. then the absolute z 
#   component is my distance to the baseplate
#
    t=len(aea)
#
#   Number of variables is 4 quaternion coefficients (put them at the end), 
#   and then an additional 'slack' variable per triangle which serves as an overhang indicator
#   the length of the array of triangle areas is thus used
    n = t + 4
#
#   The number of constraints is 1 per triangle; one to set the overhanging indicator; 
#   plus one more to constrain the quaternion coefficients to length 1
    m = t + 1
#
#   Indicator variables, followed by 4 quaternion coefficients
    x_l = np.ones(n,dtype=float)*0.
    x_u = np.ones(n,dtype=float)
    x_k = np.ones(n,dtype=float)*1.
    x_t = "I"*t+"C"*4
    x_l[-4:] = -1.
#
    x_k[-4] = 1.
    x_k[-3] = 1.
    x_k[-2] = 1.
    x_k[-1] = 1.
#
    tmp=np.load('./glob_19.npz')
#   tmp=np.load('/home/dirk/RECIPE/knap/stl/glob_48.npz')
    x_k[:] = tmp['x_i']
#
#   Last constraint is the quaternion norm (equality)
    c_s = np.ones(m,dtype=int)
    c_s[-1] = 0
    x_d = np.ones(m,dtype=float)*1e6
#
#   Pack data returned from init into aux
    aux=[ply,nrm,aea,cen,nog,t,fln,0]
#
    return n,m,x_l,x_u,x_k,x_t,x_d,c_s,aux
#
def simu(n,m,x,aux,g):
#
    f = np.zeros((m + 1), dtype=float)
    df = [np.zeros(n,dtype=float)]
#
    q = x[-4:].copy()
    [ply,nrm,aea,cen,nog,t,fln,0]=aux
    [rrm,dndr,dndi,dndj,dndk]=orien.orien_simu(q,nrm)
#
#   x = [y,q]
#
    f[0] = np.sum(x[:t]*aea)/np.sum(aea)
#
    tmp = rrm.copy()
    tmp[:,0]=-2.*tmp[:,0]; tmp[:,1]=-2.*tmp[:,1]; tmp[:,2]= 2.*tmp[:,2]
    dfdr = tmp*dndr
    dfdi = tmp*dndi
    dfdj = tmp*dndj
    dfdk = tmp*dndk
#
#   build direction
    b = np.array([0.,0.,1.]); b=b/np.linalg.norm(b)
#   minimum allowable normal component projected unto / with respect to build direction (always negative)
    a = np.cos(135/180*np.pi)
    for j in range(t):
        df[0][j] = aea[j]/np.sum(aea)
        f[j+1] =  np.dot(b,rrm[j])/a - 1. + (1./a+1.)*x[j]
        df.append((j,j,(1./a+1.)))
        df.append((j,n-4,np.dot(b,dndr[j])/a))
        df.append((j,n-3,np.dot(b,dndi[j])/a))
        df.append((j,n-2,np.dot(b,dndj[j])/a))
        df.append((j,n-1,np.dot(b,dndk[j])/a))
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
#   over=np.where(rrm[:,2]<-1./np.sqrt(2),1.0,0)
#   fd=np.sum(over*aea)/np.sum(aea)
#   print('F: %14.7e Q: (%4.1e,%4.1e,%4.1e,%4.1e)'%(fd,q[0],q[1],q[2],q[3]))
#
    return f, df
#
def apar(n):
#
    mov=1.0*np.ones(n)
    asf=[0.7,1.1]
#
    enf='t-r'
#       
    kmx=1000
    cnv=[1e-3,1e-3,1e4,1e-6,1e-3]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_t, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[np.zeros_like(x_k)]
#   if k > 0:
#       sph = f_1 - f_k
#       sph[0]=sph[0]-np.dot(df_k[0],x_1-x_k)
#       for df in df_k[1:]: sph[df[0]+1]=sph[df[0]+1]-df[2]*(x_1[df[1]]-x_k[df[1]])
#       sph=2.*sph/np.maximum(np.linalg.norm(x_1-x_k)**2.,1e-6)
#       c_x[0]=np.maximum(c_x[0]*sph[0],1e-6)
#       for j in range(len(f_k)-1):
#           for i in range(len(x_k)):
#               c_x.append((j,i,sph[j+1]))
#   else:
#       for j in range(len(f_k)-1):
#   j=len(f_k)-2
#   c_x.append((j,len(x_k)-1,2.))
#   c_x.append((j,len(x_k)-2,2.))
#   c_x.append((j,len(x_k)-3,2.))
#   c_x.append((j,len(x_k)-4,2.))
#
    d_l= np.maximum(x_l, x_k-mov*(x_u-x_l))
    d_u= np.minimum(x_u, x_k+mov*(x_u-x_l))
#
    L=L_k
    U=U_k
#
    return c_x,mov,L,U,d_l,d_u
#
