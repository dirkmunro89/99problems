#
import numpy as np
#
from prob.grnd.pNtop_swei1 import init
from prob.grnd.pNtop_swei1 import simu
from cmls import t2pl as caml
from subs.t2dual import t2d as subs
#
def apar(n):
#   
    mov=1e-1*np.ones(n,dtype=float)
    asf=[0.7,1.2]
#
    enf='none'
#
    kmx=5000
#   inf(dX), Euc(dX), dF/F, inf(KKT), viol.
    cnv=[1e-3,1e0,1e0,1e0,1e0]
#
    return mov, asf, enf, kmx, cnv
#
