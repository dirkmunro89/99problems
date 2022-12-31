#
import numpy as np
#
def init(g):
#
    n = 8; m = 16
    x_k = 400 * np.ones((n), dtype=float)
    x_l = 100 * np.ones_like(x_k)
    x_u = 1e8 * np.ones_like(x_k)
#
    c_s=np.ones(m)
#
    aux=[]
#
    return n,m,x_l,x_u,x_k,c_s,aux
#
def apar(n):
#   
    mov=2.*np.ones(n)
    asf=[3/4,4/3]
#  
    enf='none'
#     
    kmx=15
    cnv=[1e-6,1e-6,1e-6,1e-6,1e-6]
#       
    return mov, asf, enf, kmx, cnv
#
def caml(k, x_k, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    if k<=1:
        L = 0e0 #x_k-(x_u-x_l)
        U = 5e0*x_k #x_k+(x_u-x_l)
    else:
        L=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k - asf[0]*(x_1 - L_k), x_k - asf[1]*(x_1 - L_k))
        U=np.where((x_k-x_1)*(x_1-x_2) < 0e0, x_k + asf[0]*(U_k - x_1), x_k + asf[1]*(U_k - x_1))
#
    L = np.minimum(np.maximum(-50.*x_k,L),0.4*x_k)
    U = np.minimum(np.maximum(2.5*x_k,U),50*x_k)
#
    d_l = np.maximum(np.maximum(x_k/mov, 1.01*L),x_l)
    d_u = np.minimum(np.minimum(mov*x_k, 0.99*U),x_u)
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
    n_elm=8; n_nds=9; n_dim=3
#
    x_p=x
#
    n_cor=np.zeros((n_nds,n_dim),dtype=float)
    n_cor[0]=np.array([-250.,-250.,0.]); n_cor[1]=np.array([-250.,250.,0.])
    n_cor[2]=np.array([250.,250.,0.]); n_cor[3]=np.array([250.,-250.,0.])
    n_cor[4]=np.array([0.,0.,375.]); n_cor[5]=np.array([-375.,0.,0.])
    n_cor[6]=np.array([0.,375.,0.]); n_cor[7]=np.array([375.,0.,0.])
    n_cor[8]=np.array([0.,-375.,0.])
#
    e_con={}
    e_con[0]=[0,4]; e_con[1]=[1,4]; e_con[2]=[2,4]; e_con[3]=[3,4]
    e_con[4]=[5,4]; e_con[5]=[6,4]; e_con[6]=[7,4]; e_con[7]=[8,4]
#
    A = x; E = np.ones_like(A)
#
    F = np.zeros((n_dim*n_nds,1))
    fx=40e3; fy=20e3; fz=200e3
    F[12]=fx; F[13]=fy; F[14]=fz
#
    B = np.array(range(n_dim*n_nds))
    B = np.delete(B,[12, 13, 14])
#
    S={}; L={}; D={}; T={}
    K=np.zeros((n_dim*n_nds,n_dim*n_nds),dtype=float)
    sig=np.zeros(n_elm,dtype=float)
    eps=np.zeros(n_elm,dtype=float)
#
#   assembly
    for e in range(n_elm):
#
#       length
        L[e] = np.linalg.norm(n_cor[e_con[e][1]]- n_cor[e_con[e][0]])
        cxx = (n_cor[e_con[e][0]][0] - n_cor[e_con[e][1]][0])/L[e]
        cyy = (n_cor[e_con[e][0]][1] - n_cor[e_con[e][1]][1])/L[e]
        czz = (n_cor[e_con[e][0]][2] - n_cor[e_con[e][1]][2])/L[e]
#
#       d_loc=[u_2_x, u_2_y, u_2_z, u_1_x, u_1_y, u_1_z]^T
#
#       D matrix (gen def given disp. in local system)
        D[e] = np.array([[1., 0., 0., -1., 0., 0.]])
#
#       partial transformation matrix (global to local); only for truss element
        T[e] = np.zeros((6,6))
        T[e][0][0]=cxx
        T[e][0][1]=cyy
        T[e][0][2]=czz
        T[e][3][3]=cxx
        T[e][3][4]=cyy
        T[e][3][5]=czz
#
#       material law (spring stiffness, given definition of gen def)
        S[e] = A[e]*E[e]/L[e]
#
#       element stiffness matrix, transformed to global
        tmp=np.dot(D[e],T[e])
        k = S[e]*np.dot(tmp.transpose(), tmp)
#
#       assemble
        ids = np.array([n_dim*e_con[e][0],n_dim*e_con[e][0]+1,n_dim*e_con[e][0]+2, \
            n_dim*e_con[e][1],n_dim*e_con[e][1]+1,n_dim*e_con[e][1]+2])
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
    U = np.zeros(n_dim*n_nds)
    free=list(set(range(n_dim*n_nds))-set(B))
    U[free]=U0
#
#   post process
    for e in range(n_elm):
        ids = np.array([n_dim*e_con[e][0],n_dim*e_con[e][0]+1,n_dim*e_con[e][0]+2, \
            n_dim*e_con[e][1],n_dim*e_con[e][1]+1,n_dim*e_con[e][1]+2])
        dglo=U[ids]
        dloc= np.dot(T[e],dglo)
#       gen def
        eps[e]= np.dot(D[e],dloc)
#       gen stress (WHICH IS A LOAD!!!)
        sig[e]=  S[e]*eps[e]/A[e]
#
    f=np.zeros((1+m),dtype=np.float64)
    for e in range(n_elm):
        f[0]=f[0]+x_p[e]*L[e]/128211.
        f[e+1] = sig[e] -100e0
        f[e+1+8] = -sig[e] -100e0
#
    return f
#
