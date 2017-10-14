# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------
import numpy as np
from sklearn.neighbors import NearestNeighbors
from IPython import embed

def EvaluateKLcoverage(traj,opt,erg):
	xmin=opt.xlb[0]
	xmax=opt.xub[0]
	ymin=opt.xlb[1]
	ymax=opt.xub[1]
	Lx=xmax-xmin
	Ly=ymax-ymin
	if opt.sensorMode=='smallFootPrint':

		nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(opt.kdOBJ.data)
		distances, ind = nbrs.kneighbors(np.matrix(traj[:,0:2]))	
		traj_stat=opt.traj_stat +np.reshape(np.bincount(ind.T[-1],minlength=opt.traj_stat.size),(Lx,Ly),'F')
		traj_stat=traj_stat/traj_stat.sum()
		opt.traj_stat_temp=traj_stat
		

	if opt.sensorMode=='largeFootPrint':
		nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(opt.kdOBJ.data)
		distances, ind = nbrs.kneighbors(np.matrix(traj[:,0:2]))	
		opt.traj_stat_temp=np.reshape(np.bincount(ind.T[-1],minlength=opt.traj_stat.size),(Lx,Ly),'F')
		
		traj=np.matrix(traj).T	
		#gaussian processes function


		traj_stat=opt.traj_stat+np.reshape(fmu,(Lx,Ly),'F')
		traj_stat=traj_stat/traj_stat.sum()
	temp=np.reshape(opt.utility,(opt.utility.size,1),'F')	
	KL_dist=(np.multiply(temp,np.log(temp+np.finfo(float).eps))-np.multiply(temp,np.log(np.reshape(traj_stat,(traj_stat.size,1),'F')+np.finfo(float).eps))).sum()

	return KL_dist,opt,erg