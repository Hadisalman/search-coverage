import numpy as np
from math import sqrt

def accumulate_CK(opt,erg,xs):
	xmin=opt.xlb[0]
	xmax=opt.xub[0]
	ymin=opt.xlb[1]
	ymax=opt.xub[1]
	Lx=xmax-xmin
	Ly=ymax-ymin

	Nkx=erg.muk.shape[0]
	Nky=erg.muk.shape[1]

	Nsteps=xs.shape[1]
	
	temp=np.matrix((np.append([1],sqrt(0.5)*np.matrix(np.ones((1,Nky-1))))))
	HK=np.multiply(sqrt(Lx*Ly),temp.T*temp)
	KX=np.matrix(range(0,Nkx)).T*np.ones((1,Nky))
	KY=np.ones((Nkx,1))*np.matrix(range(0,Nky))
	Ck=np.zeros((erg.Nkx,erg.Nky))
	for it in range(0,Nsteps):
		xrel=xs[0,it]-xmin
		yrel=xs[1,it]-ymin

		Ck=Ck+np.divide((np.multiply(np.cos((KX)*(pi)*(xrel/Lx)),np.cos(KY*pi*yrel/Ly)))*opt.dt,HK.T)
	erg.Ck=erg.Ck+Ck	
	return erg