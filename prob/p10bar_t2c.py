#
#   todo
#
import numpy as np
#
def init(g):
#
    n = 10; m = 20
    x_k = 5. * np.ones((n), dtype=float)
    x_l = .1 * np.ones_like(x_k)
    x_u = 15. * np.ones_like(x_k)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,aux
#
def apar(n):
#   
    mov=1.0*np.ones(n)
    asf=[0.7,1.1]
# 
    enf='none'
#      
    kmx=8
    cnv=[1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, df_k, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    L=np.zeros_like(x_k)
    U=np.zeros_like(x_k)
#
#   T2C
    c_x=np.where(df_k<0,2e0*np.absolute(df_k)/x_k,0e0)
#
    c_x[0]=1e-6*np.ones_like(x_k)
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def simu(n,m,x,aux,g):
#
    f=np.zeros((1+m),dtype=float)
    df = np.zeros((m + 1, n), dtype=float)
#
    f=func(n,m,x)
#
    dx=1e-4
    for i in range(n):
        x0 = x[i]
        x[i] += dx
        fd = func(n,m,x)
        x[i] = x0
        df[:, i] = (fd - f) / dx
#
    return f, df
#
def func(n,m,x):
#
    f = np.zeros((m + 1), dtype=float)
#
    L=360
#
    n_elm=10; n_nds=6; n_dim=2
#
#   nodal coordinates
    n_cor=np.zeros((n_nds,n_dim),dtype=np.float64)
    n_cor[0]=np.array([2*L,L])
    n_cor[1]=np.array([2*L,0.])
    n_cor[2]=np.array([L,L])
    n_cor[3]=np.array([L,0.])
    n_cor[4]=np.array([0.,L])
    n_cor[5]=np.array([0.,0.])
#
#   element connectivities
    e_con={}
    e_con[0]=[4,2]
    e_con[1]=[2,0]
    e_con[2]=[5,3]
    e_con[3]=[3,1]
    e_con[4]=[2,3]
    e_con[5]=[0,1]
    e_con[6]=[4,3]
    e_con[7]=[2,5]
    e_con[8]=[2,1]
    e_con[9]=[0,3]
#
#   area and younds modulis (element props)
    A = x; E = np.ones_like(A)
#
#   global load vector
    F = np.zeros((2*n_nds,1))
    F[3]=-100
    F[7]=-100
#
#   fixed displacement boundary conditions
    B = np.zeros(4,dtype=int)
    B[0]=8
    B[1]=9
    B[2]=10
    B[3]=11
#
    S={}; L={}; D={}; T={}
    K=np.zeros((2*n_nds,2*n_nds),dtype=np.float64)
    sig=np.zeros(n_elm,dtype=np.float64)
    eps=np.zeros(n_elm,dtype=np.float64)
#
#   assembly
    for e in range(n_elm):
#
#       length
        L[e] = np.linalg.norm(n_cor[e_con[e][1]]- n_cor[e_con[e][0]])
#
#       d_loc=[u_2_x, u_2_y, u_1_x, u_1_y]^T
#
#       D matrix (gen def given disp. in local system)
        D[e] = np.array([[1., 0., -1., 0.]])#/L[e]
#
#       transformation matrix (global to local)
        cos = (n_cor[e_con[e][1]][0]- n_cor[e_con[e][0]][0])/L[e]
        sin = (n_cor[e_con[e][1]][1]- n_cor[e_con[e][0]][1])/L[e]
        T[e] = np.array([[cos, sin, 0, 0],[-sin, cos, 0, 0],[0,0,cos,sin],[0.,0.,-sin,cos]])
#
#       material law (spring stiffness, given definition of gen def)
        S[e] = A[e]*E[e]/L[e]
#
#       element stiffness matrix, transformed to global
        tmp=np.dot(D[e],T[e])
        k = S[e]*np.dot(tmp.transpose(), tmp)
#
#       assemble
        ids = np.array([2*e_con[e][1],2*e_con[e][1]+1,2*e_con[e][0],2*e_con[e][0]+1])
        K[ids[:,np.newaxis],ids]+=k
#
#   remove fixed (at zero) dofs
    tmp = np.delete(K,list(B),0)
    K0 = np.delete(tmp,list(B),1)
    F0 = np.delete(F,list(B))
#
#   solve at free dofs for loads
    U0 = np.dot(np.linalg.inv(K0),F0)
#
#   pack into complete dof array
    U = np.zeros(2*n_nds)
    free=list(set(range(2*n_nds))-set(B))
    U[free]=U0
#
#   post process
    for e in range(n_elm):
        ids = np.array([2*e_con[e][1],2*e_con[e][1]+1,2*e_con[e][0],2*e_con[e][0]+1])
        dglo=U[ids]
        dloc= np.dot(T[e],dglo)
#       gen def
        eps[e]= np.dot(D[e],dloc)
#       gen stress (WHICH IS A LOAD!!!)
        sig[e]=  S[e]*eps[e]
#
#   make response functions
    for e in range(n_elm):
        f[0]=f[0]+A[e]*L[e]*0.1
        f[e+1] = sig[e]/A[e] - 25e0
        f[e+1+n_elm] = -sig[e]/A[e] - 25e0
        if e == 8:
            f[e+1] = sig[e]/A[e] - 75e0
            f[e+1+n_elm] = -sig[e]/A[e] - 75e0
#
    return f
#
