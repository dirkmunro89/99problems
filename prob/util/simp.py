#
import numpy as np
#
def xPena(muc,penal,xPhys):
    return muc*xPhys + (1-muc)*xPhys**penal
def dxPena(muc,penal,xPhys):
    return muc + penal*(1-muc)*xPhys**(penal-1.)
#

