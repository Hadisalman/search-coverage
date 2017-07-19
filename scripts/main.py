#!/usr/bin/env python

import numpy as np
from collections import namedtuple
import time
from GenerateUtilityMap import GenerateUtilityMap
from GetFourierCoeff import GetFourierCoeff
from SMC_Update import SMC_Update
from IPython import embed
from math import pi
#from Calculate_Ergodicity import Calculate_Ergodicity
#from initialization import initialization

def initialization():
	
	pose=namedtuple('positions',['x','y','theta'])
	erg=namedtuple('ergodicity',['Nkx','Nky','mu','muk','HK','KX','KY','s'])
	opt=namedtuple('options',['Lx','Ly','nagents','Nsteps','dt','vlb','vub','wlb','wub','xmin','xmax','ymin','ymax'])
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
	erg.s=1.5
	erg.Nkx=50
	erg.Nky=50

	erg.KX=(np.matrix(range(0,erg.Nkx)).T)*np.ones((1,erg.Nky))
	erg.KY=np.ones((erg.Nkx,1))*np.matrix(range(0,erg.Nky))
	erg.LK=np.divide(1,np.power(1+(np.square(erg.KX)+np.square(erg.KY)),erg.s))

	#simulation parameters
	opt.Nsteps=5000
	opt.dt=0.1
	
	return pose,erg,opt


if __name__ == "__main__":

	pose,erg,opt=initialization()
	#addnoise=0

	X,Y,informationMap=GenerateUtilityMap(opt.xmin,opt.xmax,opt.ymin,opt.ymax)

	erg.mu=np.reshape(informationMap,X.shape)
	erg.mu=erg.mu/(erg.mu).sum()

	erg.HK,erg.muk=GetFourierCoeff(erg.mu,erg.Nkx,erg.Nky,opt.Lx,opt.Ly,X,Y)

	dt=opt.dt
	Nsteps=opt.Nsteps
	Ck=np.matrix(np.zeros((erg.Nkx,erg.Nky)))
	Ergodicity_Metric=np.matrix(np.zeros((Nsteps,1)))
	
	
	
	t=time.time()
	for it in range(0,Nsteps):
		ti=(it+1)*dt
		pose,Ck=SMC_Update(pose,opt,erg,ti,Ck)

		ck=Ck/opt.nagents/ti
		Ergodicity_Metric=(np.matrix(np.multiply(erg.LK,np.square(ck-erg.muk)))).sum()
		
	
	print("It took:")
	print(time.time()-t)	





