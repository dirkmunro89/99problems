#
import numpy as np
#
class stub:
#   instantiate with copies of data (dont return pointers)
    def __init__(self,k,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x):
        self._k = k
        self._x_k = x_k.copy()
        self._x_d = x_d.copy()
        self._mov = mov
        self._d_l = d_l.copy()
        self._d_u = d_u.copy()
        self._g_k = g_k.copy()
        self._dg_k = dg_k.copy()
        self._L_k = L_k.copy()
        self._U_k = U_k.copy()
        self._c_x = c_x.copy()
#   change the move limit and bounds of the current instantiation
    def chg(self,fct,x_l,x_u):
        mov=self._mov*fct
        self._d_l = np.maximum(self._x_k-mov*(x_u-x_l),x_l)
        self._d_u = np.minimum(self._x_k+mov*(x_u-x_l),x_u)
        self._mov=mov
        return mov
#   retrieve a copy of the current instantiationn (dont return pointers)
    def get(self):
        k = self._k
        x_k = self._x_k.copy()
        x_d = self._x_d.copy()
        d_l = self._d_l.copy()
        d_u = self._d_u.copy()
        g_k = self._g_k.copy()
        dg_k = self._dg_k.copy()
        L_k = self._L_k.copy()
        U_k = self._U_k.copy()
        c_x = self._c_x.copy()
        return k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x
#
