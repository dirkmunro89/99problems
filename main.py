#
import sys
import time
import random
import numpy as np
import logging as log
from joblib import Parallel, delayed
#
from prbs import mods
from stub import Stub
from enfc import Enfc
#
def loop(init,apar,simu,caml,subs,g):
#
    t0=time.time()
    if g > 0: print("Start %d .."%g); log = open('hist_%d.log'%g,'w')
    else: log = open('history.log','w')
#
    [n,m,x_l,x_u,x_k,aux]=init(g)
#
    [mov,asf,enf,kmx,cnv]=apar(n)
#
    if g >= 0:
        for i in range(n): x_k[i]=random.uniform(x_l[i]+.0*(x_u[i]-x_l[i]),x_u[i]-.0*(x_u[i]-x_l[i]))
    if g == -9:
        for i in range(n): x_k[i]=x_l[i]+0.5*(x_u[i]-x_l[i])
#
    x_k[:] = np.maximum(np.minimum(x_k,x_u),x_l)
#
    mov0=mov.copy(); k=0; h=[]; d_xi=1; d_xe=1; x_d=np.ones(m,dtype=float)*1e6
    x_i=x_k.copy(); x_0=x_k.copy(); x_1=x_k.copy(); x_2=x_k.copy()
    L_k=np.zeros_like(x_k); U_k=np.zeros_like(x_k)
    L=np.zeros_like(x_k); U=np.zeros_like(x_k); c_x=np.zeros((m,n))
#
    if g == -9: 
        fdck(simu,n,m,x_k,aux,0)
        return
#
    log.write(('%4s%3s%10s%12s%9s%9s%13s%8s%11s%11s\n')%\
        ('k', 's', 'Obj', 'Vio', 'Bou', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))#,flush=True)
    if not g > 0:
        print(('%4s%3s%10s%12s%9s%9s%13s%8s%11s%11s')%\
            ('k', 's', 'Obj', 'Vio', 'Bou', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))#,flush=True)
#
    enfc = Enfc()
#
    inn=0; tot=0; ts1=0.; ts0=0.; tf1=0.; tf0=0.; to=0.; to0=0.; to1=0.
    while k<kmx:
#
        ts0=time.time()
        if k > 0: f_1 = f_k.copy(); df_1 = df_k.copy()
        [f_k,df_k] = simu(n,m,x_k,aux,0); tot=tot+1
        v_k=max(f_k[1:])
        ts1=time.time(); ts=ts1-ts0
#
        tf0=time.time()
        cont=True; test=' '
        if enf == 't-r':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                if cont:
                    mov=mov*1.1
                    enfc.par_add(f_k[0],v_k,k)
                else:
                    mov=stub.set_mov(0.5,x_l,x_u)
                    [s_k,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
        elif enf == 'c-a':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.con_pas(f_1,f_k,q_k,c_x)
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
        tf1=time.time(); tf=tf1-tf0
#
        ti=time.time()
        if not str(test).strip(): itr=str(cont)[0]
        else: itr=str(test)[0]
        h.append(list(f_k)); bdd=np.count_nonzero(x_k-x_l<1e-3)/n+np.count_nonzero(x_u-x_k<1e-3)/n
        if k>0: d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)#/float(n)
        log.write('%4d%3s%14.3e%9.0e%7.2f%11.1e%11.1e%11.1e%11.1e%11.1e\n'%\
            (k, itr, f_k[0], v_k, bdd, d_xi, d_xe,ts,to,ti-to0)); log.flush()
        if not g > 0:
            print('%4d%3s%2d%14.3e%9.0e%7.2f%11.1e%11.1e%11.1e%11.1e%11.1e'%\
                (k, itr, inn, f_k[0], v_k, bdd, d_xi, d_xe, ts,to,ti-to0))#,flush=True)
#
        if k>0 and cont : 
            if d_xi<cnv[0] or d_xe<cnv[1]: 
                log.write('Termination on Convergence criteria\n')
                if not g > 0: print('Termination on Convergence criteria')
                break
        if k>0 and np.amax(mov)<cnv[0]:
            log.write('Enforced Termination; excessively reduced trust-region\n')
            if not g > 0: print('Enforced Termination; excessively reduced trust-region')
            break
        if k>0 and enf=='c-a' and inn>10:
            log.write('Enforced Termination; excessive conservatism\n')
            if not g > 0: print('Enforced Termination; excessive conservatism')
            break
#
        to0=time.time()
        if cont:
            [c_x,m_k,L,U,d_l,d_u]=caml(k,x_k,df_k,x_1,x_2,L_k,U_k,x_l,x_u,asf,mov)
            mov[:]=m_k; L_k[:]=L; U_k[:]=U; inn=0
            if enf == 't-r' or enf == 'c-a': stub=Stub(k,x_k,x_d,mov,d_l,d_u,f_k,df_k,L_k,U_k,c_x)
        else: inn=inn+1; k=k-1
#
        x_0[:]=x_k 
        [x,d,q_k] = subs(n,m,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x)
        x_k[:]=x; x_d[:]=d
        to1=time.time(); to=to1-to0; ti=time.time()
#
        if cont:
            x_2[:]=x_1; x_1[:]=x_0 
#
        k=k+1
#
    [f_k,_] = simu(n,m,x_k,aux,g)
    if g > 0: 
        np.savez_compressed('glob_%d.npz'%g, x_i=x_i, x_k=x_k, f_k=f_k); print("... exit %d"%g)
        log.close()
        return k,tot,h[-1][0],max(h[-1][1:])
    else: 
        log.write('Total number of system evaluations: %d\n'%tot); log.close()
        print('Total number of system evaluations: %d\n'%tot)
        enfc.par_plt(); enfc.cnv_plt(h)
        return h
#
def fdck(simu,n,m,x_k,aux,g):
#
#   f=np.zeros((1+m),dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    [f0,df0] = simu(n,m,x_k,aux,g)
#
    dx=1e-4
#
    print("")
    print("Error in computed derivatives with respect to finite differences")
    print("")
#
    err=-1e8
    print('%10s %10s'%("Variables", "Responses"))
    for i in range(0,n,int(np.ceil(n/100))):
        x0 = x_k[i]
        x_k[i] += dx
        [fd,_] = simu(n,m,x_k,aux,g)
        x_k[i] = x0
        df[:, i] = (fd - f0) / dx
        print("%10d "%i,end="")
        for j in range(m+1):
            print("%7.0e "%(df[j,i]-df0[j,i]),end="")
        print("")
        err=max(err,np.amax(np.absolute(df[:,i]-df0[:,i])))
    print("")
    print("Maximum absolute error: %7.0e"%err)
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
    gmx=0
    pus=0
    fdc=0
#
    if fdc:         #check finite differences
        h=loop(init,apar,simu,caml,subs,-9)
        sys.exit()
    if gmx == 0:    #standard run
        h=loop(init,apar,simu,caml,subs,-1)
    elif gmx == 1:  #one multi-start (debugging)
        h=loop(init,apar,simu,caml,subs,0)
    elif gmx>1:     #multi-starts
        fopt=1e8; gopt=0; ktot=0; stot=0; fnd=0
        res = Parallel(n_jobs=pus)(delayed(loop)(init,apar,simu,caml,subs,g) for g in range(1,gmx+1))
        glog = open('global.log','w'); g=1
        for r in res:
            (k_s,t_s,f_s,v_s)=r
            ktot=ktot+k_s
            stot=stot+t_s
            if f_s < fopt and v_s<1e-3: fopt=f_s; gopt=g; nopt='T'
            else: nopt='F'
            glog.write('%4d%3s%10d%10d%14.3e%9.0e\n'%(g, nopt, k_s, t_s, f_s, v_s))#,flush=True)
            print('%4d%3s%10d%10d%14.3e%9.0e'%(g, nopt, k_s, t_s, f_s, v_s))#,flush=True)
            g=g+1
        print("See solution %d"%gopt)
        glog.write("See solution %d\n"%gopt)
        for r in res:
            (k_s,t_s,f_s,v_s)=r
            if (f_s-fopt)/fopt < 1e-2 and v_s<1e-3: fnd=fnd+1
        print("Found %d times with %d iterations spent in total, "%(fnd,ktot)+\
            "consisting of %d system evaluations\n"%(stot))
        glog.write("Found %d times with %d iterations spent in total, "%(fnd,ktot) +\
            "consisting of %d system evaluations\n"%(stot))
        glog.close()
#
