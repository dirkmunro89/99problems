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
from wrapper import init, wrapper1, wrapper2

########################################################################################################
### MAIN FUNCTION                                                                                    ###
########################################################################################################

def main():
    # Logger
    path = os.path.dirname(os.path.realpath(__file__))
    file = os.path.join(path, "MMA_BEAM2.log")
    logger = setup_logger(file)
    logger.info("Started\n")
    # Set numpy print options
    np.set_printoptions(precision=4, formatter={'float': '{: 0.4f}'.format})
    #################################################################################################
    [n,m,x_l,x_u,x_k,aux]=init(-1); m = 1
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
    move = 0.2
    c = 1000*eeem
    d = eeem.copy()
    a0 = 1
    a = zerom.copy()
    outeriter = 0
    maxoutit = 10000
    kkttol = 1e-12 # not used
    simuc = 0
    # Calculate function values and gradients of the objective and constraints functions
    if outeriter == 0:
        f0val,df0dx,fval,dfdx = wrapper2(n,xval,aux); simuc = simuc + 1
        innerit = 0
        outvector1 = np.array([outeriter, innerit, f0val, fval])
        outvector2 = xval.flatten()
        # Log
        logger.info("outvector1 = {}".format(outvector1))
        logger.info("outvector2 = {}\n".format(outvector2))
    # The iterations starts
    kktnorm = kkttol+10
    outit = 0
    while (outit < maxoutit):
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
        # Log
        logger.info("outvector1 = {}".format(outvector1))
        logger.info("outvector2 = {}".format(outvector2))
        logger.info("kktnorm    = {}\n".format(kktnorm))

        np.savetxt('iter_%d.dat'%outeriter,xval)
        bdd=np.count_nonzero(xval-xmin<1e-3)/n+np.count_nonzero(xmax-xval<1e-3)/n
        bdd=bdd-np.count_nonzero(xmax-xmin<1e-3)/n
        ubd=np.where(xval-xmin>1e-3,np.where(xmax-xval>1e-3,1,0),0)
        mykktnorm=np.linalg.norm((df0dx + np.dot(dfdx.T,lam))*ubd,np.inf)
        mykktnorm2=np.linalg.norm((df0dx + np.dot(dfdx.T,lam))*ubd)
        f0norm=abs((f0val-f0valold)/f0val)
        print('|KKT| %e and ||KKT|| %e and |dX| %e and ||dX|| %e and BW %1.2f\n'\
            %(mykktnorm,mykktnorm2,np.linalg.norm(xval-xold1,np.inf),\
            np.linalg.norm(xval-xold1),bdd))
        print('|dF/F| %e \n'%f0norm)
        print('|viol| %e \n'%max(fval,0e0) )
        if mykktnorm < 1e-4 and np.linalg.norm(xval-xold1)<1e-1 and np.linalg.norm(xval-xold1,np.inf) < 1e-1 and f0norm < 1e-4: break

    # Final log
    logger.info(" Finished")
    #
    print('Iteration counter ', outeriter)
    print('Simulation counter ', simuc)

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

# Beam function
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
