import numpy as np
from numpy import matrix
from numpy import linalg
from math import pi,exp,sqrt

def GenerateUtilityMap(opt,addnoise):
	
	xdel=1
	ydel=1
	xRange=matrix(range(opt.xmin,opt.xmax+1-xdel,xdel))
	yRange=matrix(range(opt.ymin,opt.xmax+1-ydel,ydel))
	X=np.zeros((opt.xmax,opt.xmax))
	Y=np.zeros((opt.xmax,opt.xmax))
	
	for i in range(0,opt.ymax):
		X[i,:]=xRange
	for i in range(0,opt.xmax):
		Y[:,i]=yRange	


	m1=matrix([100.0,200.0])
	s1=matrix(150.0*np.identity(2))
	m2=matrix([220.0,200.0])
	s2=matrix(800.0*np.identity(2))
	m3=matrix([120.0,120.0])
	s3=matrix(600.0*np.identity(2))

	G1=matrix(np.zeros((X.size,1)))
	G2=matrix(np.zeros((X.size,1)))
	G3=matrix(np.zeros((X.size,1)))
	
	temp1=matrix(np.reshape(X,X.size))
	temp2=matrix(np.reshape(Y,Y.size))
	x=matrix(np.zeros((X.size,2)))
	x[:,0]=temp1.T
	x[:,1]=temp2.T

	for i in range(0,X.size):
		G1[i]=(1/(pi*sqrt(linalg.det(s1))))*(exp(-0.5*(x[i,:]-m1)*((linalg.inv(s1))*((x[i,:]-m1).T))))
		G2[i]=(1/(pi*sqrt(linalg.det(s2))))*(exp(-0.5*(x[i,:]-m2)*((linalg.inv(s2))*((x[i,:]-m2).T))))
		G3[i]=(1/(pi*sqrt(linalg.det(s3))))*(exp(-0.5*(x[i,:]-m3)*((linalg.inv(s3))*((x[i,:]-m3).T))))

	G=matrix(G1+3*G2+G3)
	G[G<0]=0	
	G=G/np.max(G)

	return X,Y,G