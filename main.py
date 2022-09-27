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
from enfc import Enfc
#
def loop(init,apar,simu,caml,subs):
#
    [n,m,x_l,x_u,x_k,aux]=init()
#
    x_k[:] = np.maximum(np.minimum(x_k,x_u),x_l)
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
    enfc = Enfc()
#
    while k<kmx:
#
        if k > 0: g_1 = g_k.copy(); dg_1 = dg_k.copy()
        [g_k,dg_k] = simu(n,m,x_k,aux)
        v_k=max(g_k[1:])
#
        cont=True; test=' '
        if enf == 't-r':
            if k == 0: enfc.par_add(g_k[0],v_k,k)
            else:
                cont=enfc.par_pas(g_1[0],g_k[0],v_k,dq)
                if cont:
                    mov=mov*1.1
                    enfc.par_add(g_k[0],v_k,k)
                else:
                    mov=stub.set_mov(0.7,x_l,x_u)
                    [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=stub.get()
        elif enf == 'c-a':
            if k == 0: enfc.par_add(g_k[0],v_k,k)
            else:
                cont=enfc.con_pas(g_1,g_k,q_k)
                if cont:
                    test=enfc.par_pas(g_1[0],g_k[0],v_k,dq)
                    if test: enfc.par_add(g_k[0],v_k,k)
                else:
                    [s_k,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x]=stub.get()
                    c_x[:]=stub.set_crv(2.,g_k,q_k)
        else:
            if k == 0: enfc.par_add(g_k[0],v_k,k)
            else: 
                test=enfc.par_pas(g_1[0],g_k[0],v_k,dq)
                if test: enfc.par_add(g_k[0],v_k,k)
        v_k=max(g_k[1:])
#
        if not str(test).strip(): itr=str(cont)[0]
        else: itr=str(test)[0]
        h.append(list(g_k)); bdd=np.count_nonzero(x_k-x_l<1e-3)/n+np.count_nonzero(x_u-x_k<1e-3)/n
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
        log.info('%3d%3s%14.3e%9.0e%7.2f%11.1e%11.1e'%\
            (k, itr, g_k[0], v_k, bdd, d_xi, d_xe))#,flush=True)
#
        if cont and k>0: 
            if d_xi<cnv[0] or d_xe<cnv[1]: break
        if mov<1e-8: break
#
        if cont:
            [c_x,L,U,d_l,d_u] = caml(k,x_k,dg_k,x_1,x_2,L_k,U_k,x_l,x_u,asf,mov)
            L_k[:]=L; U_k[:]=U
            if enf == 't-r' or enf == 'c-a': stub=Stub(k,x_k,x_d,mov,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
#
        x_0[:]=x_k 
        [x,d,dq,q_k] = subs(n,m,x_k,x_d,d_l,d_u,g_k,dg_k,L_k,U_k,c_x)
        x_k[:]=x; x_d[:]=d
#
        if cont:
            x_2[:]=x_1; x_1[:]=x_0 
#
        k=k+1
#
    enfc.par_plt()
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
