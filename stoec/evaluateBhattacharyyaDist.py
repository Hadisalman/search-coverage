import numpy as np
from math import log

def evaluateBhattacharyyaDist(traj_stat,mu):
	BC=np.sqrt(np.sqrt(np.multiply(np.reshape(traj_stat,(traj_stat.size,1),order='F'),np.reshape(mu,(mu.size,1),order='F')))).sum()

	return -log(BC)