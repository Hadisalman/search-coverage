import numpy as np
from math import pi,sqrt
			 
def GetFourierCoeff(mu,Nkx,Nky,Lx,Ly,X,Y):
	"""Function that calculates the Fourier Coefficients

	Args:
		1.mu(matrix): Information map
		2.Nkx(integer):
		3.Nky(integer):
		4.Lx(integer):
		5.Ly(integer):
		6.X(matrix):
		7.Y(matrix):

	Returns:
		1.HK(matrix):
		2.muk(matrix):

	"""	
	muk=np.matrix(np.zeros((Nkx,Nky)))

	temp=np.matrix((np.append([1],sqrt(0.5)*np.matrix(np.ones((1,Nky-1))))))
	HK=np.multiply(sqrt(Lx*Ly),temp.T*temp)

	for kx in range(0,Nkx):
		for ky in range(0,Nky):
			muk[kx,ky]=(np.matrix(np.multiply(np.multiply(mu,np.cos(kx*pi*X/Lx)),np.cos(ky*pi*Y/Ly)))).sum()/HK[ky,kx]


	return	HK,muk	