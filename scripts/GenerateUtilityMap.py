import numpy as np
from numpy import linalg
from math import pi,exp,sqrt
from scipy.stats import multivariate_normal
		
def GenerateUtilityMap(xmin,xmax,ymin,ymax,addnoise=0):
"""Function that generates the utility map
	
	Args:
		xmin(integer):
		xmax(integer):
		ymin(integer):
		ymax(integer):

	Returns:
		X(matrix):
		Y(matrix):
		G(matrix):the information map 

"""			
	xdel=1
	ydel=1
	xRange=np.matrix(range(xmin,xmax+1-xdel,xdel))
	yRange=np.matrix(range(ymin,xmax+1-ydel,ydel))
	X=np.zeros((xmax,xmax))
	Y=np.zeros((xmax,xmax))
	
	for i in range(0,ymax):
		X[i,:]=xRange
	for i in range(0,xmax):
		Y[:,i]=yRange	


	m1=[100.0,200.0]
	s1=150.0*np.identity(2)
	m2=[220.0,200.0]
	s2=800.0*np.identity(2)
	m3=[120.0,120.0]
	s3=600.0*np.identity(2)

	G1=np.matrix(np.zeros((X.size,1)))
	G2=np.matrix(np.zeros((X.size,1)))
	G3=np.matrix(np.zeros((X.size,1)))
	
	temp1=np.matrix(np.reshape(X,X.size))
	temp2=np.matrix(np.reshape(Y,Y.size))
	x=np.matrix(np.zeros((X.size,2)))
	x[:,0]=temp1.T
	x[:,1]=temp2.T
	
	mvn1 = multivariate_normal(m1,s1)
	G1 = mvn1.pdf(x)
	mvn2 = multivariate_normal(m2,s2)
	G2 = mvn2.pdf(x)
	mvn3 = multivariate_normal(m3,s3)
	G3 = mvn3.pdf(x)

	G=np.matrix(G1+3*G2+G3)
	G[G<0]=0	
	G=G/np.max(G)

	return X,Y,G