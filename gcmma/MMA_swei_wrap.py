#
import numpy as np
from prob.pNtop_swei_enf1d import init
from prob.pNtop_swei_enf1d import simu
#
def wrapper2(n,xval,aux):
    [f_k,df_k] = simu(n,1,np.reshape(xval,(n,)),aux,0) 
    f0val=f_k[0]; fval=f_k[1]
    df0dx=np.reshape(df_k[0],(n,1)); dfdx=np.array([df_k[1]])
    return f0val,df0dx,fval,dfdx
#
def wrapper1(n,xmma,aux):
    [f_k,df_k] = simu(n,1,np.reshape(xmma,(n,)),aux,0)
    f0valnew=np.array([[f_k[0]]]); fvalnew=np.array([[f_k[1]]])
    return f0valnew,fvalnew
