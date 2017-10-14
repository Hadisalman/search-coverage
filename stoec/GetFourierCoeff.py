# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------

import numpy as np
from math import pi,sqrt
from IPython import embed

def GetFourierCoeff(mu,Nkx,Nky,Lx,Ly,X,Y):


	muk=np.matrix(np.zeros((Nkx,Nky)))

	temp=np.matrix((np.append([1],sqrt(0.5)*np.matrix(np.ones((1,Nky-1))))))
	HK=np.multiply(sqrt(Lx*Ly),temp.T*temp)
	mu=np.reshape(mu,X.shape,'F')
	
	for kx in range(0,Nkx):
		for ky in range(0,Nky):
			muk[kx,ky]=(np.matrix(np.multiply(np.multiply(mu,np.cos(kx*pi*X/Lx)),np.cos(ky*pi*Y/Ly)))).sum()/HK[ky,kx]
	
	return	muk