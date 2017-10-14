# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------
from sklearn.neighbors import NearestNeighbors
import numpy as np
from math import sqrt,pi
from Calculate_Ergodicity import Calculate_Ergodicity
from IPython import embed
def EvaluateErgodicity(xs,opt,erg):

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
		#Updating Fourier Coefficients of Coverage Distribution
		xrel=xs[0,it]-xmin
		yrel=xs[1,it]-ymin

		Ck=Ck+np.divide((np.multiply(np.cos((KX)*(pi)*(xrel/Lx)),np.cos(KY*pi*yrel/Ly)))*opt.dt,HK.T)
	erg.Ck=erg.Ck+Ck
	Ergodicity_Metric=Calculate_Ergodicity(Ck/(Nsteps*opt.currentStage*opt.dt*opt.nagents),erg.muk,opt.nagents)
	
	#used in traj_cost for obstacle avoidance
	traj=xs
	nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(opt.kdOBJ.data)
	distances, ind = nbrs.kneighbors(np.matrix(traj[:,0:2]))	

	traj_stat=opt.traj_stat +np.reshape(np.bincount(ind.T[-1],minlength=opt.traj_stat.size),(Lx,Ly),'F')
	traj_stat=traj_stat/traj_stat.sum()
	opt.traj_stat_temp=traj_stat

	return Ergodicity_Metric,opt,erg