#
import sys
#
import logging as log
lvl=log.INFO; fmt='%(message)s'; hdl = [log.FileHandler('history.log',mode='w'),log.StreamHandler()]
log.basicConfig(level=lvl,format=fmt,handlers=hdl)
#
import numpy as np
from prbs import mods
from stub import Stub
#
import matplotlib.pyplot as plt
#
class Prto:
#
    def __init__(self):
        self.frnt = []
        self.sgma = 1e-1
        self.gama = 1e-5
        self.beta =  1.-self.gama
#
    def pas(self,f_1,f_k,v_k,dq):
        beta=self.beta
        gama=self.gama
        sgma=self.sgma
        if f_k+gama*v_k<=min([p[0] for p in self.frnt]) or v_k<=beta*min([p[1] for p in self.frnt]):
            df = f_1 - f_k
            if df < sgma*dq and dq > 0: # approx has descended / predicted a descent, 
            #but change in actual obj is less than 10% of predicted decrease, 
            #then something not cool, we restore, and try again
                return False
            else: # either predicted an ascent (h-type), or descent in df is sufficient (f-type)
                return True
        else: return False
#
    def add(self,f_k,v_k):
        frnt=self.frnt
        frnt.append((f_k,v_k)); frnt=sorted(frnt,key=lambda x:-x[0]); tmp=[]
        for p in frnt:
            if f_k >= p[0] or v_k >= p[1]: tmp.append(p)
        self.frnt=tmp
#
    def plt(self,f_k,v_k):
        plt.clf()
        x=[itm[1] for itm in self.frnt]
        y=[itm[0] for itm in self.frnt]
        plt.plot(x, y, color='g', linewidth=0.5)
        plt.xlabel("Vio.")
        plt.ylabel("Obj.")
        plt.scatter(x,y,marker='+',s=15,linewidths=0.5,color='b')
        plt.scatter(v_k,f_k,marker='o',s=15,linewidths=0.5,color='r',zorder=10)
        plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
        plt.grid()
        plt.savefig('prto.eps')
        plt.close()
#
def loop(init,apar,simu,caml,subs):
#
    [n,m,x_l,x_u,x_k,aux]=init()
#
    [mov,asf,enf,kmx,cnv]=apar()
#
    k=0; h=[]; d_xi=1; d_xe=1; x_d=np.ones(m,dtype=float)*1e6
    x_0=x_k.copy(); x_1=x_k.copy(); x_2=x_k.copy()
    L_k=np.zeros_like(x_k); U_k=np.zeros_like(x_k)
    L=np.zeros_like(x_k); U=np.zeros_like(x_k)
#
    log.info(('%3s%3s%14s%9s%7s%11s%11s')%\
        ('k', 's', 'Obj', 'Vio', 'Bou', '|dX|', '||dX||'))#,flush=True)
#
    prto = Prto()
#
    while k<kmx:
#
        if k > 0: g_1 = g_k.copy(); dg_1 = dg_k.copy()
        [g_k,dg_k] = simu(n,m,x_k,aux)
        v_k=max(g_k[1:])
#
        cont=True
        if enf:
            if k == 0: prto.add(g_k[0],v_k)
            else:
                cont=prto.pas(g_1[0],g_k[0],v_k,dq)
                if cont:
                    mov=mov*1.1
                    prto.add(g_k[0],v_k)
                else:
                    mov=stub.chg(0.7,x_l,x_u)
                    [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=stub.get()
            v_k=max(g_k[1:])
#
        h.append(list(g_k)); bdd=np.count_nonzero(x_k-x_l<1e-6)+np.count_nonzero(x_u-x_k<1e-6)
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
        log.info('%3d%3s%14.3e%9.0e%7i%11.1e%11.1e'%\
            (k, str(cont)[0], g_k[0], v_k, bdd, d_xi, d_xe))#,flush=True)
        if cont and k>0: 
            if d_xi<cnv[0] or d_xe<cnv[1]: break
#
        x_0[:]=x_k 
#
        if cont:
            [c_x,L,U,d_l,d_u] = caml(k, x_k, dg_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov)
            L_k[:]=L; U_k[:]=U
            if enf: stub=Stub(k,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
#
        [x,d,dq] = subs(n,m,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
        x_k[:]=x; x_d[:]=d
#
        if cont:
            x_2[:]=x_1; x_1[:]=x_0 
#
        k=k+1
#
    if enf: prto.plt(g_k[0],v_k)
#
    return h
#
def main(prob):
#
    log.info(prob)
    [init,apar,simu,caml,subs]=mods(prob)
    h=loop(init,apar,simu,caml,subs)
    return h
#
if __name__ == "__main__":
#
    if len(sys.argv)>1: prob = sys.argv[1]
    else: prob='user'
    [init,apar,simu,caml,subs]=mods(prob)
    h=loop(init,apar,simu,caml,subs)
#
