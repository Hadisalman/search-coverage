import numpy as np
from sklearn.neighbors import NearestNeighbors

def timeAverageStatisticsDistribution(traj,opt):
	nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(opt.kdOBJ)
	distances, ind = nbrs.kneighbors(np.matrix(traj[0:2,:]).T)	
	traj_stat=np.bincount(ind,minlength=opt.utility.size)
	traj_stat=traj_stat/traj_stat.sum()

	return traj_stat