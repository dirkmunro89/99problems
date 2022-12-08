function [xPhyp]=xPena(muc,penal,xPhys);
    xPhyp = muc*xPhys + (1-muc)*xPhys.^penal;
