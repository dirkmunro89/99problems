#
import numpy as np
#
def t2pl(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    tmp1=np.where(df_k[0] > 0, 3., -1.)
    tmp2=np.where(x_1*x_k>0, np.log(x_1/x_k), 0.)
    tmp3=np.where(df_k[0]*df_1[0]>0, np.log(df_1[0]/df_k[0]), 0.)
    a=np.where(tmp2*tmp3>0,1.+tmp3/tmp2,tmp1)
    c_x[0]=(a-1e0)*(df_k[0])/np.maximum(x_k,1e-6)
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t2dl(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
#   tmp1 = df_k[0] + x_d[0]*df_k[1]
#   tmp2= np.where(x_k<1e-6, 0, 1)
#   tmp3= np.where(x_k>1-1e-6, 0, 1)
    L = L_k*0.75 + np.power(df_k[0]+x_d[0]*df_k[1],2.)*0.25
#
    c_x[0] = np.sqrt(L)
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t1(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t2rl(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    c_x[0]=-2e0*(df_k[0])/np.maximum(x_k,1e-6)
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t2el(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    if k > 0:
        c_x[0]= (df_1[0] - df_k[0])/np.maximum(np.absolute(x_1-x_k),1e-6)*np.sign(x_1-x_k)
    else:
        c_x[0]= np.zeros_like(c_x[0])
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t2sl(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    if k > 0:
        sph = 2.*(f_1 - f_k - np.dot(df_k,(x_1-x_k)))/np.maximum(np.linalg.norm(x_1-x_k)**2.,1e-6)
        c_x[0]=sph[0]
    else:
        c_x[0]= np.zeros_like(c_x[0])
#
    c_x[0]=np.maximum(c_x[0],1e-6)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def t2ml(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    if k<=1:
        L = x_k-0.5*(x_u-x_l)
        U = x_k+0.5*(x_u-x_l)
    else:
        osc=(x_k-x_1)*(x_1-x_2)
        fac=np.ones_like(x_k)
        fac[np.where(osc>0)] = asf[1]
        fac[np.where(osc<0)] = asf[0]
        L=x_k-fac*(x_1-L_k)
        U=x_k+fac*(U_k-x_1)
#
    L_l = x_k - 10.*(x_u-x_l)
    L_u = x_k - 0.01*(x_u-x_l)
    U_l = x_k + 0.01*(x_u-x_l)
    U_u = x_k + 10.*(x_u-x_l)
    L=np.maximum(np.minimum(L, L_u), L_l)
    U=np.minimum(np.maximum(U, U_l), U_u)
#
    c_x[0] = np.where( df_k[0] < 0., -2./(x_k - L), 2./(U - x_k))*df_k[0]
    c_x[0] = np.maximum(c_x[0],1e-6)
#
    d_l = np.maximum(x_l, x_k - mov*(x_u-x_l))
    d_u = np.minimum(x_u, x_k + mov*(x_u-x_l))
#
    return c_x,mov,L,U,d_l,d_u
#
def mml(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
    c_x[0] = np.maximum(c_x[0],1e-6)
#
    if k>2:
        mov=np.where((x_k-x_1)*(x_1-x_2) <= 0., mov*asf[0], mov*asf[1])
#
    mov=np.minimum(np.maximum(mov,0.5e-3),.1)
#
    L=L_k
    U=U_k
#
    d_l = np.maximum(x_k-mov*(x_u-x_l),x_l)
    d_u = np.minimum(x_k+mov*(x_u-x_l),x_u)
#
    return c_x,mov,L,U,d_l,d_u
#
def mma(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    if k<=1:
        L = x_k-0.5*(x_u-x_l)
        U = x_k+0.5*(x_u-x_l)
    else:
        osc=(x_k-x_1)*(x_1-x_2)
        fac=np.ones_like(x_k)
        fac[np.where(osc>0)] = asf[1]
        fac[np.where(osc<0)] = asf[0]
        L=x_k-fac*(x_1-L_k)
        U=x_k+fac*(U_k-x_1)
#
    L_l = x_k - 10.*(x_u-x_l)
    L_u = x_k - 0.01*(x_u-x_l)
    U_l = x_k + 0.01*(x_u-x_l)
    U_u = x_k + 10.*(x_u-x_l)
    L=np.maximum(np.minimum(L, L_u), L_l)
    U=np.minimum(np.maximum(U, U_l), U_u)
#
    d_l = np.maximum(x_l, np.maximum(L + 0.1*(x_k-L), x_k - mov*(x_u-x_l)))
    d_u = np.minimum(x_u, np.minimum(U - 0.1*(U-x_k), x_k + mov*(x_u-x_l)))
#
    return c_x,mov,L,U,d_l,d_u
#
def gcm(k, x_k, x_d, f_k, df_k, f_1, df_1, x_1, x_2, L_k, U_k, x_l, x_u, asf, mov):
#
    c_x=np.zeros_like(df_k)
#
    for j in range(len(f_k)):
        c_x[j] = c_x[j] + np.sum(0.1/len(x_k)*np.absolute(df_k[j])*(x_u-x_l))
        c_x[j]=np.maximum(c_x[j],1e-6)/(x_u-x_l)
#
    if k<=1:
        L = x_k-0.5*(x_u-x_l)
        U = x_k+0.5*(x_u-x_l)
    else:
        osc=(x_k-x_1)*(x_1-x_2)
        fac=np.ones_like(x_k)
        fac[np.where(osc>0)] = asf[1]
        fac[np.where(osc<0)] = asf[0]
        L=x_k-fac*(x_1-L_k)
        U=x_k+fac*(U_k-x_1)
#
    L_l = x_k - 10.*(x_u-x_l)
    L_u = x_k - 0.01*(x_u-x_l)
    U_l = x_k + 0.01*(x_u-x_l)
    U_u = x_k + 10.*(x_u-x_l)
    L=np.maximum(np.minimum(L, L_u), L_l)
    U=np.minimum(np.maximum(U, U_l), U_u)
#
    d_l = np.maximum(x_l, np.maximum(L + 0.1*(x_k-L), x_k - mov*(x_u-x_l)))
    d_u = np.minimum(x_u, np.minimum(U - 0.1*(U-x_k), x_k + mov*(x_u-x_l)))
#
    return c_x,mov,L,U,d_l,d_u
#
