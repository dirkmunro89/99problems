#
import numpy as np
#
from prob.grnd.pNtop_comp2 import init
from prob.grnd.pNtop_comp2 import simu
from cmls import mma as caml
from subs.mma02dual import mma02 as subs
#
def apar(n):
#   
    mov=2e-1*np.ones(n,dtype=float)
    asf=[0.7,1.2]
#
    enf='none'
#
    kmx=5000
    cnv=[1e-3,1e0,1e0,1e0,1e0]
#
    return mov, asf, enf, kmx, cnv
#
