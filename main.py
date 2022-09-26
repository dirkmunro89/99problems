#
import sys
#
import logging as log
lvl=log.INFO; fmt='%(message)s'; hdl = [log.FileHandler('history.log',mode='w'),log.StreamHandler()]
log.basicConfig(level=lvl,format=fmt,handlers=hdl)
#
import numpy as np
from prbs import mods
from stub import stub
#
import matplotlib.pyplot as plt
#
def plot_pareto_frontier(Ys, Xs, L, Yo, Xo, maxX=True, maxY=True):
#
#   Xs = Xs/(max(Xs)-min(Xs))
#   Ys = Ys/(max(Ys)-min(Ys))
#
    '''Pareto frontier selection process'''
    sorted_list = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxY)
    pareto_front = [sorted_list[0]]
    for pair in sorted_list[1:]:
        if maxY:
            if pair[1] >= pareto_front[-1][1]:
                pareto_front.append(pair)
        else:
            if pair[1] <= pareto_front[-1][1]:
                pareto_front.append(pair)
 
    '''Plotting process'''
    pf_X = [pair[0] for pair in pareto_front]
    pf_Y = [pair[1] for pair in pareto_front]
    plt.plot(pf_X, pf_Y, color='g', linewidth=0.5)
    plt.xlabel("Vio.")
    plt.ylabel("Obj.")
    plt.scatter(Xs,Ys,marker='+',s=15,linewidths=0.5,color='b')
    plt.scatter(Xo,Yo,marker='o',s=15,linewidths=0.5,color='r',zorder=10)
    padx= .05*(max(pf_X) - min(pf_X))
    pady= .05*(max(pf_Y) - min(pf_Y))
    plt.xlim(min(pf_X)-padx,max(pf_X)+padx)
    plt.ylim(min(pf_Y)-pady,max(pf_Y)+pady)

#   for i in range(len(L)):
#       plt.annotate(str(L[i]), (Xs[i],Ys[i]))

    plt.show()
#
class Prto:
#
    def __init__(self):
        self.frnt = []
#
    def pas(self,f_1,f_k,v_k,dq):
        if f_k < min([itm[0] for itm in self.frnt]) or v_k < min([itm[1] for itm in self.frnt]):
            df = f_1 - f_k
            if df < 1e-1*dq and dq > 0: # approx has descended / predicted a descent, 
            #but change in actual obj is less than 10% of predicted decrease, 
            #then something not cool, we restore, and try again
                return False
            else: # either predicted an ascent (h-type), or descent in df is sufficient (f-type)
                return True
        else:
            return False
#
    def add(self,f_k,v_k):
        frnt=self.frnt
        frnt.append((f_k,v_k)); frnt=sorted(frnt,key=lambda x:-x[0]); tmp=[]
        for p in frnt:
            if f_k >= p[0] or v_k >= p[1]: tmp.append(p)
        self.frnt=tmp
#
def loop(init,apar,simu,caml,subs):
#
    [n,m,x_l,x_u,x_k,aux]=init()
#
    [mov,asf,kmx,cnv]=apar()
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
#   prto=[]
#
    rest=0
    hist=[]
#
    while k<kmx:
#
        if k > 0: g_1 = g_k.copy(); dg_1 = dg_k.copy()
        [g_k,dg_k] = simu(n,m,x_k,aux)
        v_k=max(g_k[1:])
#
        hist.append((g_k[0],v_k,k))
        hist=sorted(hist,key=lambda x:-x[0])
#
#       if acceptable to filter
        rest=0
        if k == 0: prto.add(g_k[0],v_k)#prto.append((g_k[0],v_k,k))
        else:
            if prto.pas(g_1[0],g_k[0],v_k,dq):
                prto.add(g_k[0],v_k)
                mov=mov*1.1
            else:
                mov=sub1.chg(0.7,x_l,x_u)
                [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=sub1.get()
                v_k=max(g_k[1:])
                rest=1
#           if g_k[0] < min([p[0] for p in prto]) or v_k < min([p[1] for p in prto]):
#               df = g_1[0]-g_k[0]
#               if df < 1e-1*dq and dq > 0: # approx has descended / predicted a descent, 
#               #but change in actual obj is less than 10% of predicted decrease, 
#               #then something not cool, we restore, and try again
#                   mov=sub1.chg(0.7,x_l,x_u)
#                   [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=sub1.get()
#                   v_k=max(g_k[1:])
#                   rest=2
#               else: #either predicted an ascent (h-type), or descent in df is sufficient (f-type)
#                   prto.append((g_k[0],v_k,k)); prto=sorted(prto,key=lambda x:-x[0]); tmp=[]
#                   for p in prto:
#                       if g_k[0] >= p[0] or v_k >= p[1]: tmp.append(p)
#                   prto=tmp; 
#                   mov=mov*1.1
#           else: #if not accepted by filter, restore with changed move-limit and try again
#               mov=sub1.chg(0.7,x_l,x_u)
#               [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=sub1.get()
#               v_k=max(g_k[1:])
#               rest=1
#
        h.append(list(g_k)); bdd=np.count_nonzero(x_k-x_l<1e-6)+np.count_nonzero(x_u-x_k<1e-6)
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
        log.info('%3d%3d%14.3e%9.0e%7i%11.1e%11.1e'%\
            (k, rest, g_k[0], max(g_k[1:]), bdd, d_xi, d_xe))#,flush=True)
        if rest==0 and k>0: 
            if d_xi<cnv[0] or d_xe<cnv[1]: break
#
        x_0[:]=x_k 
#
        if rest==0:
            [c_x,L,U,d_l,d_u] = caml(k, x_k, dg_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov)
            L_k[:]=L; U_k[:]=U
            sub1=stub(k,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
#
        [x,d,dq] = subs(n,m,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
        x_k[:]=x; x_d[:]=d
#
        if rest==0:
            x_2[:]=x_1; x_1[:]=x_0 #possible to move x_k=x to back and say x_1=x_0?
#
        k=k+1
#
    plt.close()
    plot_pareto_frontier(  [p[0] for p in prto] , [p[1] for p in prto] , [p[2] for p in prto], g_k[0], v_k )
    plt.savefig('prto.eps')
    plt.close()
    plot_pareto_frontier(  [p[0] for p in hist] , [p[1] for p in hist] , [p[2] for p in hist],  g_k[0], v_k, maxY=False, maxX=False )
    plt.savefig('prto2.eps')
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
