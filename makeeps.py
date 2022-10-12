#
import numpy as np
#
if __name__ == "__main__":
    nelx=2*20*3
    nely=20*3
    dat=np.load('glob_256.npz')
    x=dat['x_k']
    tmp=np.flip(x.reshape((nelx,nely)).T,0)
    tmp2=np.append(np.flip(tmp,1),tmp,axis=1)
    np.savetxt("topology.dat",tmp2,fmt='%14.7f')
