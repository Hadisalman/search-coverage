# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------

import numpy as np
from IPython import embed


def traj(z,opt,agents):
	tl=opt.tf/opt.sn
	xs=np.array(agents.xi[opt.iagent])
	t=np.array([np.arange(opt.dt,tl+opt.dt,opt.dt)])
	n_steps = t.size*opt.sn
	# z=z.flatten()
	# print(z)
	traj_sampled = np.zeros([ n_steps, xs.shape[1] ])
	

	for i in range(0,opt.sn):
			
		v=z[2*(i-1)+2]
		w=z[2*(i-1)+3]
		th=xs[-1][3]
		
		
		if opt.dim==3:
			if abs(w)<1e-10:
				traj_sampled[i*t.size:(i+1)*t.size, :] = xs + t.T*np.array([v*np.cos(th), v*np.sin(th),0,0])		

			else:
				traj_sampled[i*t.size:(i+1)*t.size, 0] = xs[0,0] + v/w * (np.sin(th + t*w) - np.sin(th))
				traj_sampled[i*t.size:(i+1)*t.size, 1] = xs[0,1] + v/w * (-np.cos(th + t*w) + np.cos(th))
				traj_sampled[i*t.size:(i+1)*t.size, 3] = xs[0,3] + t*w
		
		xs = np.array([traj_sampled[(i+1)*t.size-1, :]])

	#plotting
	# not implemented
	return traj_sampled			