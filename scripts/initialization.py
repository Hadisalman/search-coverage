import numpy as np
from numpy import matrix
from math import pi

def initialization(opt,pose):
	#Domain Bounds
	opt.xmin=0
	opt.xmax=300
	opt.ymin=0
	opt.ymax=300

	opt.Lx=opt.xmax-opt.xmin
	opt.Ly=opt.ymax-opt.ymin

	#number of agents
	opt.nagents=3

	#agents Velocity Bounds
	opt.vlb=[0.1,0.1,0.1]
	opt.vub=[5.0,5.0,5.0]

	opt.wlb=[-0.5,-0.5,-0.5]
	opt.wub=[0.5,0.5,0.5]

	#Initializing agents locations
	pose.x=[100.0,200.0,50.0]
	pose.y=[50.0,75.0,250.0]
	pose.theta=[pi/4,pi/4,pi/4]

	#ergodicity parameters
	opt.ergs=1.5
	opt.ergNkx=50
	opt.ergNky=50

	opt.ergKX=(matrix(range(0,opt.ergNkx)).T)*np.ones((1,opt.ergNky))
	opt.ergKY=np.ones((opt.ergNkx,1))*matrix(range(0,opt.ergNky))
	opt.ergLK=np.divide(1,np.power(1+(np.square(opt.ergKX)+np.square(opt.ergKY)),opt.ergs))

	#simulation parameters
	opt.Nsteps=5000
	opt.dt=0.1

	return opt,pose


