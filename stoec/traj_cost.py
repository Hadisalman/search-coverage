# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------

import numpy as np
from EvaluateErgodicity import EvaluateErgodicity
from EvaluateKLcoverage import EvaluateKLcoverage
from traj import traj
from IPython import embed

def traj_cost(z,opt,agents,erg):

	xstemp=traj(z,opt,agents)
	
	# embed()
	#1)Ergodic coverage
	if opt.algorithm=='ergodic':
		f,opt,erg=EvaluateErgodicity(xstemp[:,0:opt.dim],opt,erg)
	#2)KL divergence coverage	
	if opt.algorithm=='KL':
		f,opt,erg=EvaluateKLcoverage(xstemp[:,0:opt.dim],opt,erg)
	
	if np.multiply(opt.traj_stat_temp,opt.image).sum()>0:
		f=f+1000

	for i in range(0,2):
		
		if any(xstemp[:,i])<opt.xlb[i] or any(xstemp[:,i])>opt.xub[i]:
			f=1000
			return f

	return f	