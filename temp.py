#
import numpy as np
#
from prob.grnd.pNtop_cost1 import init
from prob.grnd.pNtop_cost1 import simu
from cmls import t2rr as caml
#from cmls import t1 as caml
from subs.t2dual import t2d as subs
#
def apar(n):
#   
    mov=2e-1*np.ones(n,dtype=float)
    asf=[0.7,1.2]
#
    enf='t-r'
#
    kmx=5000
    cnv=[1e-3,1e0,1e0,1e0,1e0]
#
    return mov, asf, enf, kmx, cnv
#
