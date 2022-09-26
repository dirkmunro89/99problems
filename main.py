#
import sys
#
debug = 0
#
from functools import partial
#
import logging as log
lvl=log.INFO; fmt='%(message)s'; hdl = [log.FileHandler('history.log',mode='w'),log.StreamHandler()]
log.basicConfig(level=lvl,format=fmt,handlers=hdl)
#
import numpy as np
from prbs import mods
#
import matplotlib.pyplot as plt
#
class stub:
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
    def chg(self,fct,x_l,x_u):
        mov=self._mov*fct
        self._d_l = np.maximum(self._x_k-mov*(x_u-x_l),x_l)
        self._d_u = np.minimum(self._x_k+mov*(x_u-x_l),x_u)
        self._mov=mov
        return mov
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
    prto=[]; rest=0
    hist=[]
#
    while k<kmx:
#
        if k > 0: g_1 = g_k.copy(); dg_1 = dg_k.copy()
        [g_k,dg_k] = simu(n,m,x_k,aux)
        v_k=max(max(g_k[1:]),-1e8)
#
        hist.append((g_k[0],v_k,k))
        hist=sorted(hist,key=lambda x:-x[0])
#
#       if acceptable to filter
        rest=0
        if k == 0: prto.append((g_k[0],v_k,k))
        else:
            print(k,g_k[0],v_k)
            if g_k[0] < min([p[0] for p in prto]) or v_k < min([p[1] for p in prto]):
                df = g_k[0] - g_1[0]
#               print('g_k',g_k[0])
#               print('g_1',g_1[0])
                if df > 1e-1*dq and dq < 0: # approx has descended / predicted a descent, 
                #and change in actual obj (hopefully also descent), 
                #is not at least 10% of predicted decrease, then something not cool, we try again
                    #htype
                    #restore
                    if 1 == debug:
                        mov=mov*0.7
                        x_k[:]=x_1
                        g_k[:]=g_1
                        dg_k[:]=dg_1
                        v_k=max(g_k[1:])
#                   print('restore1')
                    else:
                        mov=sub1.chg(0.7,x_l,x_u)
                        [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=sub1.get()
                        v_k=max(max(g_k[1:]),-1e8)
#                   print('retreiving sub from k =', s_k)
                    rest=2
                else: # update the filter and continue
                    #add to pareto front and sort such that f_j monotonically decreasing
                    prto.append((g_k[0],v_k,k)); prto=sorted(prto,key=lambda x:-x[0])
                    tmp=[]
                    #only keep pairs not dominated
                    for p in prto:
                        if g_k[0] >= p[0] or v_k >= p[1]: tmp.append(p)
                    prto=tmp
#                   print('continue')
                    mov=mov*1.1
            else: #if not accepted by filter, try again
                #restore
                if 1 == debug:
                    mov=mov*0.7
                    x_k[:]=x_1
                    g_k[:]=g_1
                    dg_k[:]=dg_1
                    v_k=max(g_k[1:])
                else:
                    mov=sub1.chg(0.7,x_l,x_u)
                    [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=sub1.get()
                    v_k=max(max(g_k[1:]),-1e8)
                    print('1 getting at ...', s_k,np.sum(x_k))
                rest=1
#
        h.append(list(g_k)); bdd=np.count_nonzero(x_k-x_l<1e-6)+np.count_nonzero(x_u-x_k<1e-6)
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
#
        log.info('%3d%3d%14.3e%9.0e%7i%11.1e%11.1e'%\
            (k, rest, g_k[0], max(g_k[1:]), bdd, d_xi, d_xe))#,flush=True)
#
        if rest==0 and k>0: 
            if d_xi<cnv[0] or d_xe<cnv[1]: break
#
        x_0[:]=x_k ## possible to remove?
#
        if 1 == debug:
            [c_x,L,U,d_l,d_u] = caml(k, x_k, dg_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov)
            L_k[:]=L; U_k[:]=U
        else:
            if rest==0:
                [c_x,L,U,d_l,d_u] = caml(k, x_k, dg_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov)
                L_k[:]=L; U_k[:]=U
                sub1=stub(k,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
                print('setting at ...', np.sum(x_k))
#
        print('...',mov,np.sum(x_k))
        [x,d,dq] = subs(n,m,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
#
        x_k[:]=x; x_d[:]=d; x_2[:]=x_1; x_1[:]=x_0 #possible to move x_k=x to back and say x_1=x_0?
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
