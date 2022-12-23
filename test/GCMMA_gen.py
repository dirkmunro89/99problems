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
Orginal work written by Krister Svanberg in Matlab. This is the python version of the code written
by Arjen Deetman. 

This script is the "beam problem" from the MMA paper of Krister Svanberg. 

    minimize 0.0624*(x(1) + x(2) + x(3) + x(4) + x(5))
    subject to 61/(x(1)^3) + 37/(x(2)^3) + 19/(x(3)^3) + 7/(x(4)^3) + 1/(x(5)^3) =< 1,
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
from MMA import gcmmasub,subsolv,kktcheck,asymp,concheck,raaupdate

# Import wrapper (Dirk)
from MMA_cwei_wrap import init,wrapper1,wrapper2

########################################################################################################
### MAIN FUNCTION                                                                                    ###
########################################################################################################

def main():
    #################################################################################################
    [n,m,x_l,x_u,x_k,c_s,aux]=init(-1); m = 1
    #################################################################################################
    epsimin = 0.0000001
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
    raa0 = 0.01
    raa = 0.01*eeem
    raa0eps = 0.000001
    raaeps = 0.000001*eeem
    outeriter = 0
    maxoutit = 8
    kkttol = 1e-12 # not used
    simuc = 0
    # Calculate function values and gradients of the objective and constraints functions
    print(('%4s%3s%8s%11s%8s%5s%12s%8s%11s%6s%9s%9s')%\
        ('k', 'l', 'Obj', 'Vio', 'Bou', 'ML', '|KKT|', '|dX|', '||dX||', 'T_s', 'T_o', 'T_t'))
    if outeriter == 0:
        f0val,df0dx,fval,dfdx = wrapper2(n,xval,aux); simuc = simuc + 1
        innerit = 0
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()
    # The iterations starts
    kktnorm = kkttol+10
    outit = 0
    simuc = 0
    hist=[]
    lam = np.ones(m)*1e6
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
        # The parameters low, upp, raa0 and raa are calculated:
        low,upp,raa0,raa= \
            asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp,raa0,raa,raa0eps,raaeps,df0dx,dfdx)
        # The MMA subproblem is solved at the point xval:
        xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp= \
            gcmmasub(m,n,iter,epsimin,xval,xmin,xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d)
        # The user should now calculate function values (no gradients) of the objective- and constraint
        # functions at the point xmma ( = the optimal solution of the subproblem).
        f0valnew,fvalnew = wrapper1(n,xmma,aux); simuc = simuc + 1
        # It is checked if the approximations are conservative:
        conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew)
        # While the approximations are non-conservative (conserv=0), repeated inner iterations are made:
        innerit = 0
        if conserv == 0:
            while conserv == 0 and innerit <= 15:
                innerit += 1
                print('%4d%3s%2d%12.3e%8.0e%6.2f%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e%9.1e'%\
                    (outeriter,'N',innerit,f0valnew,fvalnew,bdd,move,mykktnorm,\
                    np.linalg.norm(xmma-xold1,np.inf),\
                    np.linalg.norm(xmma-xold1),0,0,0))
                hist.append([f0val,fval])
                # New values on the parameters raa0 and raa are calculated:
                raa0,raa = raaupdate(xmma,xval,xmin,xmax,low,upp,f0valnew,fvalnew,f0app,fapp,raa0, \
                    raa,raa0eps,raaeps,epsimin)
                # The GCMMA subproblem is solved with these new raa0 and raa:
                xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp = gcmmasub(m,n,iter,epsimin,xval,xmin, \
                    xmax,low,upp,raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d)
                # The user should now calculate function values (no gradients) of the objective- and 
                # constraint functions at the point xmma ( = the optimal solution of the subproblem).
                f0valnew,fvalnew = wrapper1(n,xmma,aux); simuc = simuc + 1
                # It is checked if the approximations have become conservative:
                conserv = concheck(m,epsimin,f0app,f0valnew,fapp,fvalnew)
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
#
        f0norm=abs((f0val-f0valold)/f0val)
        if mykktnorm < 1e-4 and np.linalg.norm(xval-xold1)<1e-1 \
            and np.linalg.norm(xval-xold1,np.inf) < 1e-1 and f0norm < 1e-4: 
            break
    #
    print('Iteration counter ', outeriter)
    print('Simulation counter ', simuc)
    #
    return hist
########################################################################################################
### FUNCTIONS                                                                                        ###
########################################################################################################

# Setup logger
def setup_logger(logfile):
    # Create logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # Create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # Create file handler and set level to debug
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    # Add formatter to ch and fh
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # Add ch and fh to logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    # Open logfile and reset
    with open(logfile, 'w'): pass
    # Return logger
    return logger

# Beam function 1
def beam1(xval):
    nx = 5
    eeen = np.ones((nx,1))
    c1 = 0.0624
    c2 = 1
    aaa = np.array([[61.0, 37.0, 19.0, 7.0, 1.0]]).T
    xval2 = xval*xval
    xval3 = xval2*xval
    xinv3 = eeen/xval3
    f0val = c1*np.dot(eeen.T,xval)
    fval = np.dot(aaa.T,xinv3)-c2
    return f0val,fval

# Beam function 2
def beam2(xval):
    nx = 5
    eeen = np.ones((nx,1))
    c1 = 0.0624
    c2 = 1
    aaa = np.array([[61.0, 37.0, 19.0, 7.0, 1.0]]).T
    xval2 = xval*xval
    xval3 = xval2*xval
    xval4 = xval2*xval2
    xinv3 = eeen/xval3
    xinv4 = eeen/xval4
    f0val = c1*np.dot(eeen.T,xval).item()
    df0dx = c1*eeen
    fval = np.dot(aaa.T,xinv3).item()-c2
    dfdx = -3*(aaa*xinv4).T
    return f0val,df0dx,fval,dfdx


########################################################################################################
### RUN MAIN FUNCTION                                                                                ###
########################################################################################################

# Run main function / program
if __name__ == "__main__":
    main()
