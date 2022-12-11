#
import numpy as np
#
class Stub:
#   instantiate with copies of data (dont return pointers)
    def __init__(self,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x):
        self._x_k = x_k.copy()
        self._x_d = x_d.copy()
        self._mov = mov.copy()
        self._d_l = d_l.copy()
        self._d_u = d_u.copy()
        self._g_k = g_k.copy()
        self._dg_k = dg_k.copy()
        self._L_k = L_k.copy()
        self._U_k = U_k.copy()
        self._c_x = c_x.copy()
#   change the move limit and bounds of the current instantiation
    def set_mov(self,fct,x_l,x_u):
        mov=self._mov.copy()*fct
        self._d_l = np.maximum(self._x_k-mov*(x_u-x_l),x_l)
        self._d_u = np.minimum(self._x_k+mov*(x_u-x_l),x_u)
        self._mov=mov.copy()
        return mov
#   set new curvatures for old problem, and return them 
#   (because cont is false, caml update is not done, so curvatures are taken from loop)
    def set_crv(self,fct,g_k,q_k):
        c_x=self._c_x.copy()
        c_x=c_x*np.where(q_k + 1e-6 < g_k, fct , 1.)[:, np.newaxis]
        self._c_x = c_x.copy()
        return c_x
    def set_rho(self,c_x,f_k,q_k,x_k,x_0,L_k,U_k,x_u,x_l):
        dee=np.maximum(np.sum(((U_k-L_k)*(x_k-x_0)**2.)/(U_k-x_k)/(x_k-L_k)/(x_u-x_l)),1e-12)
        for j in range(len(f_k)):
            if f_k[j]>q_k[j]+0.5e-7:
                dlt=1./dee*(f_k[j]-q_k[j])/(x_u-x_l)
                c_x[j]=np.minimum(1.1*(c_x[j]+dlt),10.*c_x[j])
        self._c_x = c_x.copy()
        return c_x
#   retrieve a copy of the current instantiationn (dont return pointers)
    def get(self):
        x_k = self._x_k.copy()
        x_d = self._x_d.copy()
        d_l = self._d_l.copy()
        d_u = self._d_u.copy()
        g_k = self._g_k.copy()
        dg_k = self._dg_k.copy()
        L_k = self._L_k.copy()
        U_k = self._U_k.copy()
        c_x = self._c_x.copy()
        return x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x
#
