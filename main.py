#
import sys
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
def plot_pareto_frontier(Xs, Ys, maxX=True, maxY=True):
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
    plt.scatter(Xs,Ys)
    pf_X = [pair[0] for pair in pareto_front]
    pf_Y = [pair[1] for pair in pareto_front]
    plt.plot(pf_X, pf_Y)
    plt.xlabel("Objective 1")
    plt.ylabel("Objective 2")
    plt.show()
#
def loop(init,apar,simu,caml,subs):
#
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
        [g_k,dg_k] = simu(n,m,x_k,aux); v_k=max(g_k[1:])
#
        hist.append((g_k[0],v_k))
#
#       if acceptable to filter
        rest=0
        if k == 0: prto.append((g_k[0],v_k))
        else:
            if g_k[0] < min([p[0] for p in prto]) or v_k < min([p[1] for p in prto]):
                df = g_k[0] - g_1[0]
                if df > 1e-1*dq and dq < 0: # approx has descended / predicted a descent, 
                #and change in actual obj (hopefully also descent), 
                #is not at least 10% of predicted decrease, then something not cool, we try again
                    #htype
                    #restore
                    mov=mov*0.7
                    x_k[:]=x_1
                    g_k[:]=g_1
                    dg_k[:]=dg_1
                    v_k=max(g_k[1:])
#                   print('restore1')
                    rest=2
                else: # update the filter and continue
                    #add to pareto front and sort such that f_j monotonically decreasing
                    prto.append((g_k[0],v_k)); prto=sorted(prto,key=lambda x:-x[0])
                    tmp=[]
                    #only keep pairs not dominated
                    for p in prto:
                        if g_k[0] >= p[0] or v_k >= p[1]: tmp.append(p)
                    prto=tmp
#                   print('continue')
                    mov=mov*1.1
            else: #if not accepted by filter, try again
                #restore
                mov=mov*0.7
                x_k[:]=x_1
                g_k[:]=g_1
                dg_k[:]=dg_1
                v_k=max(g_k[1:])
#               print('restore2')
                rest=1
#
        h.append(list(g_k)); bdd=np.count_nonzero(x_k-x_l<1e-6)+np.count_nonzero(x_u-x_k<1e-6)
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
#
        log.info('%3d%3d%14.3e%9.0e%7i%11.1e%11.1e'%\
            (k, rest, g_k[0], v_k, bdd, d_xi, d_xe))#,flush=True)
#
        if rest==0 and k>0: 
            if d_xi<cnv[0] or d_xe<cnv[1]: break
#
        x_0[:]=x_k
#
#       to restore from here: 
#           x_k, dg_k, x_1, x_2, L_k, U_k, mov(!)
#           ###returned from caml: c_x, L_k, U_k, d_l, d_u
#           n,m,k,x_d,g_k
#
        [c_x,L,U,d_l,d_u] = caml(k, x_k, dg_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov)
        L_k[:]=L; U_k[:]=U
#
        [x,d,dq] = subs(n,m,k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
#
        x_k[:]=x; x_d[:]=d; x_2[:]=x_1; x_1[:]=x_0
#
        k=k+1
#
    print(prto)
    plt.close()
    plot_pareto_frontier(  [p[0] for p in prto] , [p[1] for p in prto]  )
    plt.savefig('prto.eps')
    plt.close()
    plot_pareto_frontier(  [p[0] for p in hist] , [p[1] for p in hist]  )
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
