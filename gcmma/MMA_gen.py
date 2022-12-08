########################################################################################################
### GCMMA-MMA-Python         															             ### 
###                                                                                                  ###
### This file is part of GCMMA-MMA-Python. GCMMA-MMA-Python is licensed under the terms of GNU       ###
### General Public License as published by the Free Software Foundation. For more information and    ###
### the LICENSE file, see <https://github.com/arjendeetman/GCMMA-MMA-Python>.                        ###
###                                                                                                  ###
### The orginal work is written by Krister Svanberg in MATLAB.                                       ###
### This is the python version of the code written Arjen Deetman.                                    ###
### version 09-11-2019                                                                               ###
########################################################################################################

"""
Orginal work is written by Krister Svanberg in Matlab. This is the python version of the code written
by Arjen Deetman. 

This script is the "beam problem" from the MMA paper of Krister Svanberg. 

    minimize 0.0624*(x(1) + x(2) + x(3) + x(4) + x(5))
    subject to 61/(x(1)^3) + 37/(x(2)^3) + 19/(x(3)^3) +  7/(x(4)^3) +  1/(x(5)^3) =< 1,
               1 =< x(j) =< 10, for j=1,..,5. 
"""

########################################################################################################
### LOADING MODULES                                                                                  ###
########################################################################################################

# Loading modules
from __future__ import division
import numpy as np
import logging
import sys
import os

# Import MMA functions
from MMA import mmasub,subsolv,kktcheck

# Import wrapper (Dirk)
from MMA_swei_wrap import init, wrapper1, wrapper2

########################################################################################################
### MAIN FUNCTION                                                                                    ###
########################################################################################################

def main():
    #################################################################################################
    [n,m,x_l,x_u,x_k,c_s,aux]=init(-1); m = 1
    #################################################################################################
    eeen = np.ones((n,1))
    eeem = np.ones((m,1))
    zeron = np.zeros((n,1))
    zerom = np.zeros((m,1))
    xval = np.reshape(x_k,(n,1))#5*eeen
    xold1 = xval.copy()
    xold2 = xval.copy()
    xmin = np.reshape(x_l,(n,1))#eeen.copy()
    xmax = np.reshape(x_u,(n,1))#10*eeen
    low = xmin.copy()
    upp = xmax.copy()
    move = 0.1
    c = 1000*eeem
    d = eeem.copy()
    a0 = 1
    a = zerom.copy()
    outeriter = 0
    maxoutit = 10
    kkttol = 1e-12 # not used
    simuc = 0
    print(('%4s%3s%8s%11s%8s%5s%12s%8s%11s%6s%9s%9s')%\
        ('k', 'l', 'Obj', 'Vio', 'Bou', 'ML', '|KKT|', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))
    # Calculate function values and gradients of the objective and constraints functions
    if outeriter == 0:
        f0val,df0dx,fval,dfdx = wrapper2(n,xval,aux); simuc = simuc + 1
        innerit = 0
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()
    # The iterations starts
    kktnorm = kkttol+10
    outit = 0
    lam = np.ones(m)*1e6
    hist=[]
    while (outit < maxoutit):
        bdd=np.count_nonzero(xval-xmin<1e-3)/n+np.count_nonzero(xmax-xval<1e-3)/n
        bdd=bdd-np.count_nonzero(xmax-xmin<1e-3)/n
        ubd=np.where(xval-xmin>1e-3,np.where(xmax-xval>1e-3,1,0),0)
        mykktnorm=np.linalg.norm((df0dx + np.dot(dfdx.T,lam))*ubd,np.inf)
        print('%4d%3s%2d%12.3e%8.0e%6.2f%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e'%\
            (outeriter,'N',innerit,f0val,fval,bdd,move,mykktnorm,np.linalg.norm(xval-xold1,np.inf),\
            np.linalg.norm(xval-xold1),0,0,0))
        hist.append([f0val,fval])
        outit += 1
        outeriter += 1
        # The MMA subproblem is solved at the point xval:
        xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp = \
            mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,move)
        # Some vectors are updated:
        xold2 = xold1.copy()
        xold1 = xval.copy()
        xval = xmma.copy()
        # Re-calculate function values and gradients of the objective and constraints functions
        f0valold=f0val
        f0val,df0dx,fval,dfdx = wrapper2(n,xval,aux); simuc = simuc + 1
        # The residual vector of the KKT conditions is calculated
        residu,kktnorm,residumax = \
            kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()

        f0norm=abs((f0val-f0valold)/f0val)
        if np.linalg.norm(xval-xold1,np.inf) < 1e-2: 
            break
    #
    print('Iteration counter ', outeriter)
    print('Simulation counter ', simuc)
    #
    return hist

########################################################################################################
### RUN MAIN FUNCTION                                                                                ###
########################################################################################################

# Run main function / program
if __name__ == "__main__":
    main()
