#
def xPena(muc,penal,xPhys):
    return xPhys/(1.+(penal-1.)*(1.-xPhys))
def dxPena(muc,penal,xPhys):
    return (1.+(penal-1.)*(1.-xPhys)-xPhys*(1.-penal))/(1.+(penal-1.)*(1.-xPhys))**2.
