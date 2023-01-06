#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
#
#
class Enfc:
#
    def __init__(self):
        self.pf = []
        self.sgma = 0e-2
        self.gama = 0e-4
        self.kapa = 0e-4
        self.beta = 1.-self.gama
#
    def con_pas(self,g_1,g_k,q_k,c_x):
        gama=self.gama
        f_k = g_k[0]; f_1 = g_1[0]
        # if feasible descent
        if f_k < f_1 and np.amax(g_k[1:]) < 0.: return True
        # if conservative
        else:
            if np.any(np.where(q_k < g_k,1,0)):
                return False
            return True
#
    def gcm_pas(self,g_k,q_k):
        # if conservative
        if np.any(np.where(q_k + 1e-6 < g_k, 1, 0)):
            return False
        return True
#
    def par_pas(self,f_1,f_k,v_k,p_k):
        sgma=self.sgma
        gama=self.gama
        kapa=self.kapa
        beta=self.beta
        pf=self.pf
        # check if acceptable to filter
        if f_k+gama*v_k<min([p[0] for p in pf]) or v_k<beta*max(min([p[1] for p in pf]),0.):
            df = f_1 - f_k # actual descent
            dq = f_1 - p_k # predicted descent
            if df < sgma*dq and dq > kapa*v_k**2.: 
                return False
            else:
                # if h-type step
                if dq <= kapa*v_k**2.:
                    tmp=[]
                    for p in pf:
                        if f_k > p[0] or v_k > p[1]: tmp.append((p[0],p[1]))
                    tmp.append((f_k,v_k))
                    self.pf=tmp
                return True
        else: return False
#
    def par_add(self,f_k,v_k,k):
        pf=self.pf
        pf.append((f_k,v_k,k))
        self.pf=pf
#
