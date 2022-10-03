#
import sys
import random
import numpy as np
import logging as log
#
from prbs import mods
from stub import Stub
from enfc import Enfc
#
def loop(init,apar,simu,caml,subs,g):
#
    if g > 0: log = open('hist_%d.log'%g,'w')
    else: log = open('history.log','w')
#
    [n,m,x_l,x_u,x_k,aux]=init(g)
#
    [mov,asf,enf,kmx,cnv]=apar()
#
    if g >= 0:
        for i in range(n): x_k[i]=random.uniform(x_l[i],x_u[i])
#
    x_k[:] = np.maximum(np.minimum(x_k,x_u),x_l)
#
    mov0=mov; k=0; h=[]; d_xi=1; d_xe=1; x_d=np.ones(m,dtype=float)*1e6
    x_i=x_k.copy(); x_0=x_k.copy(); x_1=x_k.copy(); x_2=x_k.copy()
    L_k=np.zeros_like(x_k); U_k=np.zeros_like(x_k)
    L=np.zeros_like(x_k); U=np.zeros_like(x_k); c_x=np.zeros((m,n))
#
    log.write(('%3s%3s%14s%9s%7s%11s%11s\n')%\
        ('k', 's', 'Obj', 'Vio', 'Bou', '|dX|', '||dX||'))#,flush=True)
    if not g > 0:
        print(('%3s%3s%14s%9s%7s%11s%11s')%\
            ('k', 's', 'Obj', 'Vio', 'Bou', '|dX|', '||dX||'))#,flush=True)
#
    enfc = Enfc()
#
    while k<kmx:
#
        if k > 0: f_1 = f_k.copy(); df_1 = df_k.copy()
        [f_k,df_k] = simu(n,m,x_k,aux,0)
        v_k=max(f_k[1:])
#
        cont=True; test=' '
        if enf == 't-r':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                if cont:
                    mov=mov*1.1
                    enfc.par_add(f_k[0],v_k,k)
                else:
                    mov=stub.set_mov(0.7,x_l,x_u)
                    [s_k,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
        elif enf == 'c-a':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.con_pas(f_1,f_k,q_k)
                if cont:
                    test=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                    if test: enfc.par_add(f_k[0],v_k,k)
                else:
                    [s_k,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
                    c_x[:]=stub.set_crv(2.,f_k,q_k)
        else:
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else: 
                test=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                if test: enfc.par_add(f_k[0],v_k,k)
        v_k=max(f_k[1:])
#
        if not str(test).strip(): itr=str(cont)[0]
        else: itr=str(test)[0]
        h.append(list(f_k)); bdd=np.count_nonzero(x_k-x_l<1e-3)/n+np.count_nonzero(x_u-x_k<1e-3)/n
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
        log.write('%3d%3s%14.3e%9.0e%7.2f%11.1e%11.1e\n'%\
            (k, itr, f_k[0], v_k, bdd, d_xi, d_xe)); log.flush()
        if not g > 0:
            print('%3d%3s%14.3e%9.0e%7.2f%11.1e%11.1e'%\
                (k, itr, f_k[0], v_k, bdd, d_xi, d_xe))#,flush=True)
#
        if k>0 and cont : 
            if d_xi<cnv[0] or d_xe<cnv[1]: 
                [f_k,_] = simu(n,m,x_k,aux,g)
                log.write('Termination on Convergence criteria\n')
                if not g > 0: print('Termination on Convergence criteria')
                break
        if k>0 and mov<cnv[0]:
            [f_k,_] = simu(n,m,x_k,aux,g)
            log.write('Enforced Termination; excessively reduced trust-region\n')
            if not g > 0: print('Enforced Termination; excessively reduced trust-region')
            break
        if k>0 and np.amax(c_x)>1e16:
            [f_k,_] = simu(n,m,x_k,aux,g)
            log.write('Enforced Termination; excessive conservatism\n')
            if not g > 0: print('Enforced Termination; excessive conservatism')
            break
#
        if cont:
            [c_x,L,U,d_l,d_u] = caml(k,x_k,df_k,x_1,x_2,L_k,U_k,x_l,x_u,asf,mov); L_k[:]=L; U_k[:]=U
            if enf == 't-r' or enf == 'c-a': stub=Stub(k,x_k,x_d,mov,d_l,d_u,f_k,df_k,L_k,U_k,c_x)
#
        x_0[:]=x_k 
        [x,d,q_k] = subs(n,m,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x)
        x_k[:]=x; x_d[:]=d
#
        if cont:
            x_2[:]=x_1; x_1[:]=x_0 
#
        k=k+1
#
    log.close()
#
    if g > 0: 
        np.savez_compressed('glob_%d.npz'%g, x_i=x_i, x_k=x_k, f_k=f_k)
        return h[-1][0],max(h[-1][1:])
    else: 
        enfc.par_plt(); enfc.cnv_plt(h)
        return h
#
def main(prob):
#
    [init,apar,simu,caml,subs]=mods(prob)
    h=loop(init,apar,simu,caml,subs,-2)
#
    return h
#
if __name__ == "__main__":
#
    if len(sys.argv)>1: prob = sys.argv[1]
    else: prob='user'
#
    [init,apar,simu,caml,subs]=mods(prob)
#
#   0 to do standard run
#   1 to do one standard run with a random start (test of mult start)
#   X to do X random multi-starts
#
    gmx=1
#
    if gmx == 0:    #standard run
        h=loop(init,apar,simu,caml,subs,-1)
    elif gmx == 1:  #one multi-start (debugging)
        h=loop(init,apar,simu,caml,subs,0)
    elif gmx>1:     #multi-starts
        fopt=1e8; gopt=0
        for g in range(1,gmx+1):
            [f_s,v_s]=loop(init,apar,simu,caml,subs,g)
            if f_s < fopt and v_s<1e-3: fopt=f_s; gopt=g
            print("done")
        print(g,fopt)
#
