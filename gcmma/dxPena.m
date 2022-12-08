function [dxPhyp]=dxPena(muc,penal,xPhys);
    dxPhyp = muc + penal*(1-muc)*xPhys.^(penal-1.);
