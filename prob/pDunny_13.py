#
import numpy as np
import importlib.util
from subs.t2milx import t2c as subs
#from subs.t2duel import t2d as subs
#from prob.util.orien import orien_init
#from prob.util.orien import orien_simu
#from prob.util.orien import orien_outp
#
spec=importlib.util.spec_from_file_location("orien", "/home/dirk/RECIPE/knap/orien_nuppt.py")
orien = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orien)
#
def init(g):
#
#   Run init from RECIPE module
    [ply,nrm,aea,cen,nog,fln]=orien.orien_init("/home/dirk/RECIPE/knap/stl/dunny.stl",0)
#
#   stl is read in, and moved to a new cog, such that a sphere which just touches each ref. coordinate 
#   plane, fits around it. The positioned vtp (with some data) is written to path/name.vtp; incl.
#   _shw, _sup, _int, and _sph
#
#   Number of triangles
    t=len(aea)
#
#   Number of variables is 4 quaternion coefficients (put them at the end), 
#   and then an additional integer 'slack' variable per triangle which serves as an overhang indicator
#   the length of the array of triangle areas is thus used
    n = t + 4 + 2
#
#   The number of constraints is 1 per triangle; one to set the overhanging indicator; 
#   plus one more to constrain the quaternion coefficients to length 1
    m = 3*t
#
#   Indicator variables, followed by 4 quaternion coefficients
    x_l = np.ones(n,dtype=float)*0.
    x_u = np.ones(n,dtype=float)
    x_k = np.ones(n,dtype=float)*1.
    x_k[0]=0.
    x_k[1]=1.
    x_t ="CC"+"I"*t+"C"*4
    x_l[-4:] = -1.
    x_u[-4:] = 1.
#
    x_k[-4] = 0e-6
    x_k[-3] = 0e-6
    x_k[-2] = 0e-6
    x_k[-1] = 1e6
#   x_k[-4:]=np.array([-0.82075593, -0.33162381,  0.45697361,  0.08695103])
#
    tmp=np.load('./prob/pDunny_13_g/glob_10.npz')
    x_k[:] = tmp['x_i']
#
#   Last constraint is the quaternion norm (equality)
    c_s = np.ones(m,dtype=int)
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
#   build direction
    b = np.array([0.,0.,1.]); b=b/np.linalg.norm(b)
#
    q = x[-4:].copy()
    [ply,nrm,aea,cen,nog,t,fln,k]=aux
    xog = nog - np.array([0.,0.,x[0]*nog[2]])
    [rrm,ren,dndr,dndi,dndj,dndk,dcdr,dcdi,dcdj,dcdk]=orien.orien_simu(q,nrm,cen-nog)
#
#   x = [c,h,y,q]
#
#   for testing
#   x[:t]=np.where(rrm[:,2]<-1./np.sqrt(2),1,0)
#
    fs_a = np.sum(aea)
    fs_b = np.sum(aea*(cen[:,2]))
    f[0] = np.sum(x[2:t+2]*aea)/fs_a
    f[0] = f[0] + np.sum(x[2:t+2]*aea*(np.dot(ren+xog,b)))/fs_b + x[1]#np.sum(np.dot(ren+xog,b))
    df[0][0]= np.sum(x[2:t+2]*aea*(np.dot(-np.array([0,0,nog[2]]),b)))/fs_b 
    df[0][1] = 1.
#   df[0][0]=df[0][0] + np.sum(np.dot(-np.array([0,0,nog[2]]),b))
#
#   minimum allowable normal component projected unto / with respect to build direction (always negative)
    a = np.cos(135/180*np.pi)
    for j in range(t):
#
        df[0][j+2] = aea[j]/fs_a + aea[j]*(ren[j,2]+xog[2])/fs_b
#
        f[j+1] =  np.dot(b,rrm[j])/a - 1. + (1./a+1.)*x[j+2]
        df.append((j,j+2,(1./a+1.)))
        df.append((j,n-4,np.dot(b,dndr[j])/a))
        df.append((j,n-3,np.dot(b,dndi[j])/a))
        df.append((j,n-2,np.dot(b,dndj[j])/a))
        df.append((j,n-1,np.dot(b,dndk[j])/a))
#
        f[t+j+1] = -ren[j][2]-xog[2] + 5.
        df.append((t+j,0,float(nog[2])))
        df.append((t+j,n-4,-np.dot(dcdr[j],b)))
        df.append((t+j,n-3,-np.dot(dcdi[j],b)))
        df.append((t+j,n-2,-np.dot(dcdj[j],b)))
        df.append((t+j,n-1,-np.dot(dcdk[j],b)))
#
        f[2*t+j+1] = ren[j][2]+xog[2]-2.*x[1]*nog[2]
        df.append((2*t+j,0,-float(nog[2])))
        df.append((2*t+j,1,-2.*float(nog[2])))
        df.append((2*t+j,n-4,np.dot(dcdr[j],b)))
        df.append((2*t+j,n-3,np.dot(dcdi[j],b)))
        df.append((2*t+j,n-2,np.dot(dcdj[j],b)))
        df.append((2*t+j,n-1,np.dot(dcdk[j],b)))
#
    df[0][-4]=np.sum(x[2:t+2]*aea*(np.dot(dcdr,b)))/fs_b
    df[0][-3]=np.sum(x[2:t+2]*aea*(np.dot(dcdi,b)))/fs_b
    df[0][-2]=np.sum(x[2:t+2]*aea*(np.dot(dcdj,b)))/fs_b
    df[0][-1]=np.sum(x[2:t+2]*aea*(np.dot(dcdk,b)))/fs_b
#
#   unit quaternion constraint and derivative terms
#    f[-1] = q[0]**2. + q[1]**2. + q[2]**2. + q[3]**2. - 1.0
#    df.append((m-1,n-4,2.*q[0]))
#    df.append((m-1,n-3,2.*q[1]))
#    df.append((m-1,n-2,2.*q[2]))
#    df.append((m-1,n-1,2.*q[3]))
#
    kt=aux[-1]
    kt=kt+1
    aux[-1]=kt
#
    if g == 0:
        orien.orien_outp(fln+'_pos.vtp',x[0]*nog[2],x[-4:],kt)
#
#   ovr_srf=np.where(rrm[:,2]<-1./np.sqrt(2)-1e-6,aea,0)
#   ovr_vol=np.where(rrm[:,2]<-1./np.sqrt(2)-1e-6,aea*(ren[:,2]+nog[2]), 0.)
#   scl=1.0
#   fd= scl*np.sum(ovr_srf)/np.sum(aea)
#   fd=fd+(scl)*np.sum(ovr_vol)/np.sum(aea*(cen[:,2]))
#   print('F: %14.7e Q: (%4.1e,%4.1e,%4.1e,%4.1e)'%(fd,q[0],q[1],q[2],q[3]))
#
    return f, df
#
def apar(n):
#
    mov=0.1*np.ones(n)
    asf=[0.7,1.1]
#
    enf='none'
#       
    kmx=1000
    cnv=[1e-3,1e-3,1e4,1e6,1e-3]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, x_t, f_k, df_k, f_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=[np.zeros_like(x_k)]
#
    i=0
    d_l=np.zeros_like(x_l)
    d_u=np.zeros_like(x_u)
    for c in x_t:
        if c =='C':
            d_l[i]= np.maximum(x_l[i], x_k[i]-mov[i]*(x_u[i]-x_l[i]))
            d_u[i]= np.minimum(x_u[i], x_k[i]+mov[i]*(x_u[i]-x_l[i]))
        else:
            d_l[i]= np.maximum(x_l[i], x_k[i]-(x_u[i]-x_l[i]))
            d_u[i]= np.minimum(x_u[i], x_k[i]+(x_u[i]-x_l[i]))
        i=i+1
#
    L=L_k
    U=U_k
#
    return c_x,mov,L,U,d_l,d_u
#
