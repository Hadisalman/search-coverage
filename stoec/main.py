a#!/usr/bin/env python

import numpy as np
from numpy import matlib
from collections import namedtuple 
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from scipy import spatial
from math import pi,sqrt,log
import time
import cv2
from IPython import embed
from GenerateUtilityMap import GenerateUtilityMap
from GetFourierCoeff import GetFourierCoeff
from evaluateBhattacharyyaDist import evaluateBhattacharyyaDist
from gen_traj_CE import gen_traj_CE
# from plot_initial import plot_initial

def initialize_gen_traj_CE():
	
	opt=namedtuple('options',['xmin','xmax','ymin','ymax','Lx','Ly','dim','xlb','xub','ng','sn','dt','tf','X','Y','Z','xss','utility','nagents'
					,'colors','traj_stat','agents','vlb','vub','wlb','wub','iagent','currentStage','figg','iters','stages','devBias','map'
					,'algorithm','mapWithObstacles','sensoeMode','image','z','kdOBJ','htraj'])
	
	ce=namedtuple('cross_entropy',['N','v','iter','sigma','C0','z0','z','C','lb','ub','Flag','Fig'])

	opt.xmin=0
	opt.xmax=150
	opt.ymin=0
	opt.ymax=150

	opt.Lx=opt.xmax-opt.xmin
	opt.Ly=opt.ymax-opt.ymin

	opt.dim=3
	opt.xlb=[opt.xmin,opt.ymin]
	opt.xub=[opt.xmax,opt.ymax]
	opt.ng=[opt.Lx,opt.Ly]
	opt.sn=5
	opt.dt=0.1
	opt.tf=40
	opt.vlb=0.1
	opt.vub=5.0
	opt.wlb=-0.2
	opt.wub=0.2
	opt.iters=1

	ce.N=40
	ce.v=0.8
	ce.iter=6
	ce.sigma=0
	
	ce.C0=np.matrix(np.diag(np.concatenate(np.matlib.repmat([1],opt.sn*2,1))))
	ce.z0=np.matlib.repmat(np.matrix([1,0]).T,opt.sn,1)

	opt.stages=10
	opt.devBias=0.001

	#need to check opt.map
	ce.Flag=0

	opt.z=ce.z0
	ce.z=ce.z0
	ce.C=ce.C0

	ce.lb=np.matlib.repmat(np.matrix([opt.vlb,opt.wlb]).T,opt.sn,1)
	ce.ub=np.matlib.repmat(np.matrix([opt.vub,opt.wub]).T,opt.sn,1)
	return opt,ce

if __name__ == "__main__":
		
	agents=namedtuple('agents',['xi','xps','trajFigOPTIMAL','trajFig'])
	gp_model=namedtuple('gaussian_processes',['mean','cov','lik','inf'])
	erg=namedtuple('ergodicity',['mu','Nkx','Nky','muk','Ck','HK'])
	
	opt,ce=initialize_gen_traj_CE()
	
	opt.algorithm = 'KL' #'KL --> KL_STOEC in the paper' - 'ergodic --> ergodic_STOEC in the paper'
	opt.sensorMode = 'smallFootPrint' #choose between 'smallFootPrint' and 'largeFootPrint' if using KL
	
	img=cv2.imread('obstacleMap.png',0) #read image
	#img=cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) #covert it from rgb to gray scale ,,its already gray(0 in cv2.imread for gray/1 for colored )
	img=cv2.normalize(img.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX) #normalize 
	img=cv2.resize(img,(opt.ng[0],opt.ng[1]), interpolation=cv2.INTER_LINEAR) 
	img=np.flipud(img)
	opt.image=img

	img2=cv2.imread('infoMapwithObstacles.png',1)
	img2=cv2.normalize(img2.astype('float'), None, 0.0, 1.0, cv2.NORM_MINMAX)
	img2=cv2.resize(img,(opt.ng[0],opt.ng[1]), interpolation=cv2.INTER_LINEAR) 
	opt.mapWithObstacles=img2
	# embed()
	X,Y,utility=GenerateUtilityMap(opt.xmin,opt.xmax,opt.ymin,opt.ymax)

	utility=np.matrix(np.reshape(utility,X.shape))
	temp=np.matrix(np.ones(img.shape)-img)
	utility=np.matrix(np.multiply(utility,temp),)
	
	opt.X=np.matrix(X)
	opt.Y=np.matrix(Y)
	opt.Z=np.matrix(np.zeros(X.shape))

	opt.xss=np.matrix(np.zeros((X.size,3))) 
	opt.xss[:,0]=np.reshape(X,(X.size,1),'F')
	opt.xss[:,1]=np.reshape(Y,(X.size,1),'F')
	opt.xss[:,2]=np.reshape(opt.Z,(X.size,1),'F')

	opt.utility=utility/utility.sum()
	opt.nagents=3
	
	agents.xi=np.matrix([[120,30,0,90*pi/180],[130,130,0,270*pi/180],[30,120,0,270*pi/180]])
	agents.xps=agents.xi
	
	opt.colors=['m','k','w']
	opt.kdOBJ=spatial.KDTree(opt.xss[:,0:2])
	#####initialize ergodicity###############
	erg.mu=opt.utility

	Nk=10
	erg.Nkx=Nk
	erg.Nky=Nk

	erg.muk=GetFourierCoeff(erg.mu,erg.Nkx,erg.Nky,opt.Lx,opt.Ly,X,Y)
	erg.Ck=np.matrix(np.zeros((Nk,Nk)))

	#####initialize GP###############
	sn=0.01
	ell=7
	sf=sqrt(1)
	Ncg=30
	gp_para_lik=log(sn)
	gp_para_cov=np.log(np.matrix([ell,sf]).T)
	fmu=np.zeros((opt.xss[:,1]).shape)
	fs2=np.zeros((opt.xss[:,1]).shape)
	Lx=opt.ng[0]
	Ly=opt.ng[1]

	opt.traj_stat=np.zeros((Lx,Ly))
	
	# opt,agents,ce.Fig=plot_initial(opt,agents)
	
	t=time.time()
	BhattDistance=np.zeros((opt.stages,1))
	save_traj_stat=np.zeros((150*150,opt.stages))
	euclidean_dist=np.zeros((opt.stages,1))
	
	full_trajectory = []
	for k in range(0,opt.stages):#####look if k from 0 or 1
		opt.currentStage=k
		
		for iagent in range(0,opt.nagents):
			opt.iagent=iagent

			opt,ce,xs=gen_traj_CE(opt,ce,agents,erg)
			
			#plotting
			# embed()
			xf=xs[-1,:]
			agents.xi[iagent]=xf
			# agents.xps[iagent]=np.concatenate(agents.xps[iagent],np.xs[1:,:])
			full_trajectory = full_trajectory + xs.tolist()
			
			# plt.figure(1)
			
			# plt.title('Optimal trajectory')
			
			
				
			
			# plt.plot(np.array(full_trajectory)[:,0],np.array(full_trajectory)[:,1])
			# plt.axis((opt.xmin,opt.xmax,opt.ymin,opt.ymax))


			if opt.algorithm=='ergodic':
				erg=accumulate_CK(opt,erg,xs)

			#if opt.sensorMode=='largeFootPrint' and opt.algorithm=='KL':
				#gp
				

			if opt.algorithm=='KL' and opt.sensorMode=='smallFootPrint' or opt.algorithm=='ergodic':
				
				nbrs = NearestNeighbors(n_neighbors=1, algorithm='ball_tree').fit(opt.kdOBJ.data)
				traj = xs
				distances, ind = nbrs.kneighbors(np.matrix(traj[:,0:2]))
				opt.traj_stat=opt.traj_stat +np.reshape(np.bincount(ind.T[-1],minlength=opt.traj_stat.size),(Lx,Ly),'F')
				#plotting
		# embed()
		traj_stat_normalized=opt.traj_stat/opt.traj_stat.sum()
		save_traj_stat[:,k]=np.reshape(traj_stat_normalized,(traj_stat_normalized.size,1),'F').T  
		euclidean_dist[k,0]=np.square(traj_stat_normalized-opt.utility).sum()
		BhattDistance[k,0]=evaluateBhattacharyyaDist(traj_stat_normalized,erg.mu)
		# plt.show()