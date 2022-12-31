from __future__ import division
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt ;import cvxopt.cholmod
#
def K(nelx,nely,rmin,penal,x,u,fixed):
#
	ft=1 # ft==0 -> sens, ft==1 -> dens
#
	# Max and min stiffness
	Emin=1e-9
	Emax=1.0
#
	# dofs:
	ndof = 2*(nelx+1)*(nely+1)
#
	# Allocate design variables (as array), initialize and allocate sens.
	xold=x.copy()
	xPhys=x.copy()
#
	# FE: Build the index vectors for the for coo matrix format.
	KE=lk()
	edofMat=np.zeros((nelx*nely,8),dtype=int)
	for elx in range(nelx):
		for ely in range(nely):
			el = ely+elx*nely
			n1=(nely+1)*elx+ely
			n2=(nely+1)*(elx+1)+ely
			edofMat[el,:]=np.array([2*n1+2, 2*n1+3, 2*n2+2, 2*n2+3,2*n2, 2*n2+1, 2*n1, 2*n1+1])
	# Construct the index pointers for the coo format
	iK = np.kron(edofMat,np.ones((8,1))).flatten()
	jK = np.kron(edofMat,np.ones((1,8))).flatten()    
#
	# Filter: Build (and assemble) the index+data vectors for the coo matrix format
	nfilter=int(nelx*nely*((2*(np.ceil(rmin)-1)+1)**2))
	iH = np.zeros(nfilter)
	jH = np.zeros(nfilter)
	sH = np.zeros(nfilter)
	cc=0
	for i in range(nelx):
		for j in range(nely):
			row=i*nely+j
			kk1=int(np.maximum(i-(np.ceil(rmin)-1),0))
			kk2=int(np.minimum(i+np.ceil(rmin),nelx))
			ll1=int(np.maximum(j-(np.ceil(rmin)-1),0))
			ll2=int(np.minimum(j+np.ceil(rmin),nely))
			for k in range(kk1,kk2):
				for l in range(ll1,ll2):
					col=k*nely+l
					fac=rmin-np.sqrt(((i-k)*(i-k)+(j-l)*(j-l)))
					iH[cc]=row
					jH[cc]=col
					sH[cc]=np.maximum(0.0,fac)
					cc=cc+1
	# Finalize assembly and convert to csc format
	H=coo_matrix((sH,(iH,jH)),shape=(nelx*nely,nelx*nely)).tocsc()	
	Hs=H.sum(1)

	# BC's and support
	dofs=np.arange(2*(nelx+1)*(nely+1))
	free=np.setdiff1d(dofs,fixed)

	# Solution and RHS vectors
	f=np.zeros((ndof,1))
	u=np.zeros((ndof,1))

	# Set load
	f[1,0]=-1
	# Initialize plot and plot the initial design
    
	h=[]

	loop=0
	change=1
	dc = np.ones(nely*nelx)
	ddc = np.ones(nely*nelx)
	ce = np.ones(nely*nelx)

	fdof = ndof - len(fixed)

	sK=((KE.flatten()[np.newaxis]).T*(Emin+(xPhys)**penal*(Emax-Emin))).flatten(order='F')
	K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
	# Remove constrained dofs from matrix and convert to coo
	K = deleterowcol(K,fixed,fixed)
	D = K.diagonal()
	D = cvxopt.spmatrix(D, range(fdof), range(fdof))
	K=K.tocoo()
	# Solve system 
	K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int)) ###!!!
	B = cvxopt.matrix(f[free,0])
#
	ce[:] = 1./2.*(np.dot(u[edofMat].reshape(nelx*nely,8),KE)*u[edofMat].reshape(nelx*nely,8)).sum(1)
	dc[:]=(penal*xPhys**(penal-1)*(Emax-Emin))*ce + np.ones(nelx*nely,dtype=float)
	ddc[:]=np.maximum((penal*(penal-1)*xPhys**(penal-2)*(Emax-Emin))*ce ,1e-9)
#
	dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0] 
	ddc[:] = np.asarray(H*(ddc[np.newaxis].T/Hs))[:,0] 
#
	sH = np.append(sK,ddc)
	iH = np.append(iK, range(ndof,ndof+nelx*nely))
	jH = np.append(jK, range(ndof,ndof+nelx*nely))
	H = coo_matrix((sH,(iH,jH)),shape=(ndof+nelx*nely,ndof+nelx*nely)).tocsc()
	H = deleterowcol(H,fixed,fixed).tocoo()
	H = cvxopt.spmatrix(H.data,H.row.astype(int),H.col.astype(int)) ###!!!
#H = H.diagonal()
#H = cvxopt.spmatrix(H, range(fdof+nelx*nely), range(fdof+nelx*nely))
#H=H.tocoo()
	G = cvxopt.matrix(  np.append(f[free,0],dc ) )
#
#U = cvxopt.matrix(f[free,0])
#print(K)
#cvxopt.cholmod.linsolve(K,U)
#
	return H, D, G
#
def lk():
	E=1
	nu=0.3
	k=np.array([1/2-nu/6,1/8+nu/8,-1/4-nu/12,-1/8+3*nu/8,-1/4+nu/12,-1/8-nu/8,nu/6,1/8-3*nu/8])
	KE = E/(1-nu**2)*np.array([ [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
	[k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
	[k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
	[k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
	[k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
	[k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
	[k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
	[k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]] ]);
	return (KE)
#
def deleterowcol(A, delrow, delcol):
	# Assumes that matrix is in symmetric csc form !
	m = A.shape[0]
	keep = np.delete (np.arange(0, m), delrow)
	A = A[keep, :]
	keep = np.delete (np.arange(0, m), delcol)
	A = A[:, keep]
	return A    
#
