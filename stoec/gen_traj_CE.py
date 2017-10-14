# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------

import numpy as np
import time
from cem import cem
from IPython import embed
from math import cos,sin
from traj_cost import traj_cost
from traj import traj

def gen_traj_CE(opt,ce,agents,erg):
	iagents=opt.iagent
	ce.C=0.5*ce.C+0.5*ce.C0
	z=ce.z0
	
	t=time.time()
	for i in range(1,opt.iters+1):
		opt.z=z
		ce.Flag=1
		
		z,_,_,C = cem(traj_cost,z,ce,opt,agents,erg)
		# embed()
		# print(z)
		ce.C=C
		# embed()
		xs=traj(z,opt,agents)
		# embed()
		# embed()
		#plotting
	
	
	print('Duration:')
	print(time.time()-t)	
	return opt,ce,xs
	