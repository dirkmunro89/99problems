# A 200 LINE TOPOLOGY OPTIMIZATION CODE BY NIELS AAGE AND VILLADS EGEDE JOHANSEN, JANUARY 2013
# Updated by Niels Aage February 2016
from __future__ import division
import numpy as np
import random
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from matplotlib import colors
import matplotlib.pyplot as plt
import cvxopt ;import cvxopt.cholmod
#
def xPena(muc,pen,xPhys):
	return muc*xPhys + (1-muc)*xPhys**pen
def dxPena(muc,pen,xPhys):
	return muc + pen*(1-muc)*xPhys**(pen-1.)
#
def main(nelx,nely,volfrac,pen,qen,muc,rmin,ft):
	print("Minimum compliance problem with OC")
	print("ndes: " + str(nelx) + " x " + str(nely))
	print("volfrac: " + str(volfrac) + ", rmin: " + str(rmin) + ", pen: " + str(pen))
	print("Filter method: " + ["Sensitivity based","Density based"][ft])

	mov=0.1*np.ones(nelx*nely)

	# Max and min stiffness
	Emin=1e-9
	Emax=1.0

	# dofs:
	ndof = 2*(nelx+1)*(nely+1)

	# Allcate design variables (as array), initialize and allocate sens.
	x= volfrac * np.ones(nely*nelx,dtype=float)
	xold=x.copy(); xolder=x.copy(); xoldest=x.copy()
	xPhys=x.copy()

	dc=np.zeros((nely,nelx), dtype=float)

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
	fixed=np.union1d(dofs[0:2*(nely+1):2],np.array([2*(nelx+1)*(nely+1)-1]))
	free=np.setdiff1d(dofs,fixed)

	# Solution and RHS vectors
	f=np.zeros((ndof,1))
	u=np.zeros((ndof,1))

	# Set load
	f[1,0]=-1*0.
	# Initialize plot and plot the initial design
	plt.ion() # Ensure that redrawing is possible
	fig,ax = plt.subplots()
	im = ax.imshow(-x.reshape((nelx,nely)).T, cmap='gray',\
	interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
	fig.show()
    
	gv=-9.81/(nelx*nely)
	loop=0
	change=1
	dv = np.ones(nely*nelx)
	dc = np.ones(nely*nelx)
	ce = np.ones(nely*nelx)
	op = 0e0
	while change>0.01 and loop<2000:
		loop=loop+1
		# Setup and solve FE problem
		sK=((KE.flatten()[np.newaxis]).T*(Emin+xPena(muc,pen,xPhys)*(Emax-Emin))).flatten(order='F')
		K = coo_matrix((sK,(iK,jK)),shape=(ndof,ndof)).tocsc()
		# Remove constrained dofs from matrix and convert to coo
		K = deleterowcol(K,fixed,fixed).tocoo()
        # Set self-weight load
		f_tmp=np.zeros(ndof)
		np.add.at(f_tmp,edofMat[:, 1::2].flatten(),np.kron(xPena(muc,qen,xPhys),gv*np.ones(4)/4.))
		f_apl = f.copy()
		f_apl[:,0] += f_tmp
		# Solve system 
		K = cvxopt.spmatrix(K.data,K.row.astype(int),K.col.astype(int))
		B = cvxopt.matrix(f_apl[free,0])
		cvxopt.cholmod.linsolve(K,B)
		u[free,0]=np.array(B)[:,0] 

		# Objective and sensitivity
		ce[:] = (np.dot(u[edofMat].reshape(nelx*nely,8),KE)*u[edofMat].reshape(nelx*nely,8) ).sum(1)
		obj=( (Emin+xPhys**pen*(Emax-Emin))*ce ).sum()
		dc[:]=(-dxPena(muc,pen,xPhys)*(Emax-Emin))*ce
		dc[:] += 2. * gv * u[edofMat[:,1],0] * 1./4. * dxPena(muc,qen,xPhys)
		dc[:] += 2. * gv * u[edofMat[:,3],0] * 1./4. * dxPena(muc,qen,xPhys)
		dc[:] += 2. * gv * u[edofMat[:,5],0] * 1./4. * dxPena(muc,qen,xPhys)
		dc[:] += 2. * gv * u[edofMat[:,7],0] * 1./4. * dxPena(muc,qen,xPhys)

		dv[:] = -np.ones(nely*nelx)
		# Sensitivity filtering:
		if ft==0:
			dc[:] = np.asarray((H*(x*dc))[np.newaxis].T/Hs)[:,0] / np.maximum(0.001,x)
		elif ft==1:
			dc[:] = np.asarray(H*(dc[np.newaxis].T/Hs))[:,0]
			dv[:] = np.asarray(H*(dv[np.newaxis].T/Hs))[:,0]

		op = op*0.1 + np.power(dc,2.)*0.9
		g=-np.sum(xPhys)/nelx/nely/volfrac+1.
		# Optimality criteria
		xold[:]=x
		if loop > 2:
			dx=x-xold
			xp = xp*0.9 + np.power(dx,2.)*0.1
			osc=(x-xolder)*(xolder-xoldest)
			mov=np.minimum(np.maximum(1e-3,np.where(osc<=0.,mov*0.7,mov*1.2)),1e-1)
		xp = np.power(np.ones_like(x),2.)
		(x[:],g)=oc(nelx,nely,x,volfrac,dc,dv,g,mov,op,xp)
		xoldest[:]=xolder; xolder[:]=xold

		# Filter design variables
		if ft==0:   xPhys[:]=x
		elif ft==1:	xPhys[:]=np.asarray(H*x[np.newaxis].T/Hs)[:,0]
	    
		# Compute the change by the inf. norm
		change=np.linalg.norm(x.reshape(nelx*nely,1)-xold.reshape(nelx*nely,1),np.inf)

		# Plot to screen
		im.set_array(-x.reshape((nelx,nely)).T)
		fig.canvas.draw()
		fig.canvas.flush_events()

		# Write iteration history to screen (req. Python 2.6 or newer)
		print("it.: {0} , obj.: {1:.3e} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
					loop,obj,(g+volfrac*nelx*nely)/(nelx*nely),change))

	# Make sure the plot stays and that the shell remains	
	plt.show()
	input("Press any key...")
    
#element stiffness matrix
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

def oc(nelx,nely,x,volfrac,dc,dv,g,mov,op,xp):
	l1=1e-9
	l2=1e9
	# reshape to perform vector operations
	xnew=np.zeros(nelx*nely)

	while (l2-l1)/(l1+l2)>1e-9:
		lmid=0.5*(l2+l1)
		c_k = np.sqrt(op)/np.sqrt(xp)#-2.*dc/x
		tmp=np.minimum(np.maximum((dc + lmid*dv)/c_k,-mov),mov)
		xnew[:] = np.minimum(np.maximum(0.,x - tmp),1.)
		gt=g+np.sum((dv*(xnew-x)))
		if gt>0:
			l1=lmid
		else:
			l2=lmid
	return (xnew,gt)
    
def deleterowcol(A, delrow, delcol):
	# Assumes that matrix is in symmetric csc form !
	m = A.shape[0]
	keep = np.delete (np.arange(0, m), delrow)
	A = A[keep, :]
	keep = np.delete (np.arange(0, m), delcol)
	A = A[:, keep]
	return A    

if __name__ == "__main__":
	# Default input parameters
	nelx=120
	nely=60
	volfrac=0.1
	rmin=3.3
	pen=3.
	qen=1.
	muc=1e-2
	ft=1 # ft==0 -> sens, ft==1 -> dens

	import sys
	if len(sys.argv)>1: nelx   =int(sys.argv[1])
	if len(sys.argv)>2: nely   =int(sys.argv[2])
	if len(sys.argv)>3: volfrac=float(sys.argv[3])
	if len(sys.argv)>4: rmin   =float(sys.argv[4])
	if len(sys.argv)>5: pen  =float(sys.argv[5])
	if len(sys.argv)>6: ft     =int(sys.argv[6])

	main(nelx,nely,volfrac,pen,qen,muc,rmin,ft)
