import os
import math
import numpy as np
c=0
vecs_0 = 0.
while os.path.isfile('u_vec_%d.dat'%c):
    vecs=np.loadtxt('u_vec_%d.dat'%c)
    if c==0:
        vecs_0=vecs
    num=np.zeros_like(vecs)
    i=0
    for v_0, v in zip(vecs_0,vecs):
        num[i]=np.dot(v_0,v)/np.linalg.norm(v)**2.
        i=i+1
    print(np.amin(num),np.amax(num))
    c=c+1
#
