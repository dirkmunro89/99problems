#
import sys
import time
import random
import inspect
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
    [n,m,x_l,x_u,x_k,c_s,aux]=init(g)
#
    if np.count_nonzero(c_s==0) > 0 and 't2duel' not in inspect.getfile(subs): 
        print('Error: equality constraints implemented in only t2d subsolver from t2duel.py'); sys.exit()
#
    [mov,asf,enf,kmx,cnv]=apar(n)
#
    if g >= 0:
        for i in range(n): x_k[i]=random.uniform(x_l[i]+.0*(x_u[i]-x_l[i]),x_u[i]-.0*(x_u[i]-x_l[i]))
    if g == -9:
        for i in range(n): x_k[i]=random.uniform(x_l[i]+.0*(x_u[i]-x_l[i]),x_u[i]-.0*(x_u[i]-x_l[i]))
#
    x_k[:] = np.maximum(np.minimum(x_k,x_u),x_l)
#
    scl0=1.;mov0=mov.copy(); k=0; h=[]; d_xi=1; d_xe=1; x_d=np.ones(m,dtype=float)*1e6
    x_i=x_k.copy(); x_0=x_k.copy(); x_1=x_k.copy(); x_2=x_k.copy()
    L_k=np.zeros_like(x_k); U_k=np.zeros_like(x_k)
    L=np.zeros_like(x_k); U=np.zeros_like(x_k); c_x=np.zeros((m+1,n)); f_1 = np.zeros(m+1)
#
    if g == -9: 
        fdck(simu,n,m,x_k,aux,0)
        return
#
    log.write(('%4s%3s%8s%12s%7s%5s%12s%8s%11s%6s%9s%9s\n')%\
        ('k', 'l', 'Obj', 'Vio', 'Bou', 'ML', '|KKT|', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))
    if not g > 0:
        print(('%4s%3s%8s%12s%7s%5s%12s%8s%11s%6s%9s%9s')%\
            ('k', 'l', 'Obj', 'Vio', 'Bou', 'ML', '|KKT|', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))
#
    enfc = Enfc()
#
    inn=0; tot=0; ts1=0.; ts0=0.; tf1=0.; tf0=0.; to=0.; to0=0.; to1=0.
    while k<kmx:
#
        ts0=time.time()
        if k > 0: f_1 = f_k.copy(); df_1 = df_k.copy()
        [f_k,df_k] = simu(n,m,x_k,aux,0); tot=tot+1
        if k == 0: scl0=f_k[0]
        f_k[0]=f_k[0]/scl0; df_k[0]=df_k[0]/scl0
        v_k=np.maximum(np.amax(f_k[1:]*c_s),0.)
        if 0 in c_s: v_k=np.maximum(np.amax(np.absolute(f_k[1:]*(1-c_s))),v_k)
        ts1=time.time(); ts=ts1-ts0
#
        tf0=time.time()
        cont=True; test=' '
        if enf == 't-r':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                if cont:
                    mov=mov*1.2; mov=np.minimum(mov,mov0)
                    enfc.par_add(f_k[0],v_k,k)
                else:
                    mov=stub.set_mov(0.7,x_l,x_u)
                    [x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
        elif enf == 'c-a':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.con_pas(f_1,f_k,q_k,c_x)
                if cont:
                    test=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                    if test: enfc.par_add(f_k[0],v_k,k)
                else:
                    [x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
                    c_x[:]=stub.set_crv(10.,f_k,q_k)
        elif enf == 'gcm':
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else:
                cont=enfc.gcm_pas(f_k,q_k)
                if cont:
                    test=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                    if test: enfc.par_add(f_k[0],v_k,k)
                else:
                    c_x[:]=stub.set_rho(c_x,f_k,q_k,x_k,x_0,L_k,U_k,x_u,x_l)
                    [x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x]=stub.get()
        else:
            if k == 0: enfc.par_add(f_k[0],v_k,k)
            else: 
                test=enfc.par_pas(f_1[0],f_k[0],v_k,q_k[0])
                if test: enfc.par_add(f_k[0],v_k,k)
        tf1=time.time(); tf=tf1-tf0
#
        ti=time.time()
        if not str(test).strip(): itr=str(cont)[0]
        else: itr=str(test)[0]
        tmp=f_k.copy(); tmp[0]=tmp[0]*scl0; h.append(list(tmp))
        bdd=np.count_nonzero(x_k-x_l<1e-3)/n+np.count_nonzero(x_u-x_k<1e-3)/n
        bdd=bdd-np.count_nonzero(x_u-x_l<1e-3)/n
        ubd=np.where(x_k-x_l>1e-3,np.where(x_u-x_k>1e-3,1,0),0)
        kkt=np.linalg.norm((df_k[0] + np.dot(x_d,df_k[1:]))*ubd,np.inf)
#
        v_k=np.amax(f_k[1:]*c_s)
        if 0 in c_s: v_k=np.maximum(np.amax(np.absolute(f_k[1:]*(1-c_s))),v_k)
#
        if k == 0: ti=to0
        else: 
            d_xi=np.linalg.norm(x_k-x_0,np.inf); d_xe=np.linalg.norm(x_k-x_0)
            d_f0=abs((f_k[0]-f_1[0])/f_k[0])
        log.write('%4d%3s%2d%12.3e%8.0e%6.2f%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e\n'%\
            (k,itr,inn,f_k[0]*scl0,v_k,bdd,np.amax(mov),kkt,d_xi,d_xe,ts,to,ti-to0)); log.flush()
        if not g > 0:
            print('%4d%3s%2d%12.3e%8.0e%6.2f%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e'%\
                (k,itr,inn,f_k[0]*scl0,v_k,bdd,np.amax(mov),kkt,d_xi,d_xe,ts,to,ti-to0))#,flush=True)
#
        if k>1 and cont:
            if d_xi<cnv[0] and d_xe<cnv[1] and d_f0<cnv[2] and kkt<cnv[3] and v_k<cnv[4]:
                log.write('Termination on Convergence criteria\n')
                if not g > 0: print('Termination on Convergence criteria')
                break
        if k>1 and enf=='t-r' and inn>30:
            log.write('Enforced Termination; excessively reduced trust-region\n')
            if not g > 0: print('Enforced Termination; excessively reduced trust-region')
            break
        if k>1 and enf=='c-a' and inn>30:
            log.write('Enforced Termination; excessive conservatism\n')
            if not g > 0: print('Enforced Termination; excessive conservatism')
            break
        if k>1 and enf=='gcm' and inn>30:
            log.write('Enforced Termination; excessive conservatism\n')
            if not g > 0: print('Enforced Termination; excessive conservatism')
            break
#
        to0=time.time()
        if cont:
            [c_x,m_k,L,U,d_l,d_u]=caml(k,x_k,f_k,df_k,f_1,x_1,x_2,L_k,U_k,x_l,x_u,asf,mov)
            mov[:]=m_k; L_k[:]=L; U_k[:]=U; inn=0
            if enf == 't-r' or enf == 'c-a' or enf == 'gcm': 
                stub=Stub(x_k,x_d,mov,d_l,d_u,f_k,df_k,L_k,U_k,c_x)
        else: inn=inn+1; k=k-1
#
        x_0[:]=x_k 
        [x,d,q_k] = subs(n,m,x_k,x_d,d_l,d_u,f_k,df_k,L_k,U_k,c_x,c_s)
        x_k[:]=x; x_d[:]=d
        to1=time.time(); to=to1-to0; ti=time.time()
#
        if cont: x_2[:]=x_1; x_1[:]=x_0 
#
        k=k+1
#
    f_k[0]=f_k[0]*scl0
    df_k[0]=df_k[0]*scl0
#
    np.savetxt('x_str.log',x_k)
    np.savetxt('d_str.log',x_d)
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
    df = np.zeros((m + 1, n), dtype=float)
#
    [f0,df0] = simu(n,m,x_k,aux,g)
    scl0=f0[0]
    f0[0]=f0[0]/scl0; df0[0]=df0[0]/scl0
#
    dx=1e-6
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
        fd[0]=fd[0]/scl0
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
        fopt=1e16; gopt=0; ktot=0; stot=0; fnd=0
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
        for r in res:
            (k_s,t_s,f_s,v_s)=r
            if (f_s-fopt)/fopt < 1e-2 and v_s<1e-3: fnd=fnd+1
        print("Solution %5d, f=%7.3e, found %6d time(s)."%(gopt,fopt,fnd))
        print("%12d iterations spent in total, with"%(ktot))
        print("%12d system evaluations,"%(stot))
        print("%12.0f per (probably) optimal solution."%(stot/fnd))
        glog.write("Solution %5d, f=%7.3e, found %6d time(s).\n"%(gopt,fopt,fnd))
        glog.write("%12d iterations spent in total, with\n"%(ktot))
        glog.write("%12d system evaluations,\n"%(stot))
        glog.write("%12.0f per (probably) optimal solution.\n"%(stot/fnd))
#
