import numpy as np
from numpy import matrix
from math import pi,sqrt

def GetFourierCoeff(opt,X,Y):

	mu=opt.ergmu
	Nkx=opt.ergNkx
	Nky=opt.ergNky

	opt.ergmuk=matrix(np.zeros((Nkx,Nky)))

	temp=matrix((np.append([1],sqrt(0.5)*matrix(np.ones((1,Nky-1))))))
	opt.ergHK=np.multiply(sqrt(opt.Lx*opt.Ly),temp.T*temp)

	for kx in range(0,Nkx):
		for ky in range(0,Nky):
			opt.ergmuk[kx,ky]=(matrix(np.multiply(np.multiply(mu,np.cos(kx*pi*X/opt.Lx)),np.cos(ky*pi*Y/opt.Ly)))).sum()/opt.ergHK[ky,kx]


	return	opt	