import numpy as np
from numpy import matrix
from math import pi,cos,sin
from IPython import embed

def SMC_Update(pose,opt,time,Ck):
	
	Lx=opt.Lx
	Ly=opt.Ly
	dt=opt.dt
	KX=opt.ergKX
	KY=opt.ergKY
	HK=opt.ergHK
	muk=opt.ergmuk
	LK=opt.ergLK

	for iagent in range(0,opt.nagents):
		xrel=pose.x[iagent]-opt.xmin
		yrel=pose.y[iagent]-opt.ymin
		
		Ck=Ck+np.divide((np.multiply(np.cos((KX)*(pi)*(xrel/Lx)),np.cos(KY*pi*yrel/Ly)))*dt,HK.T)
		
	
	temp=np.multiply(np.divide(LK,HK.T),Ck-opt.nagents*time*muk)	
	
	for j in range(0,opt.nagents):
		xRel=pose.x[j]-opt.xmin
		yRel=pose.y[j]-opt.ymin

		Bjx=(np.multiply(temp,np.multiply(np.multiply(-KX*pi/Lx,np.sin(KX*pi*xRel/Lx)),np.cos(KY*pi*yRel/Ly)))).sum()
		Bjy=(np.multiply(temp,np.multiply(np.multiply(-KY*pi/Ly,np.cos(KX*pi*xRel/Lx)),np.sin(KY*pi*yRel/Ly)))).sum()
		
		GammaV=Bjx*cos(pose.theta[j])+Bjy*sin(pose.theta[j])
		GammaW=-Bjx*sin(pose.theta[j])+Bjy*cos(pose.theta[j])
		
		if GammaV >= 0:
			v=opt.vlb[j]
		if GammaV < 0:
			v=opt.vub[j]
		if GammaW >= 0 :
			w=opt.wlb[j]
		if GammaW < 0 :
			w=opt.wub[j]
		
		if abs(w)< 1e-10:
			pose.x[j] = pose.x[j] + v*dt*cos(pose.theta[j])   
			pose.y[j] = pose.y[j] + v*dt*sin(pose.theta[j])
	   
		if abs(w)>=1e-10:
			
			pose.x[j] = pose.x[j] + v/w*(sin(pose.theta[j] + w*dt) - sin(pose.theta[j]))
			pose.y[j] = pose.y[j] + v/w*(cos(pose.theta[j]) - cos(pose.theta[j]+ w*dt))		
			

		pose.theta[j]=pose.theta[j]+ w*dt
		
	return 	pose,Ck	

			



