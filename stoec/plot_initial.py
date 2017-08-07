from matplotlib import pyplot
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from math import pi
from IPython import embed
from draw_path import draw_path
from traj import traj

def plot_initial(opt,agents):

	plt.figure(1)
	
	#plt.setp(axes,fontsize=20)
	ax2=plt.subplot(1,3,2)
	plt.title('Optimal trajectory')
	#ax=Axes3D(plt.gcf())
	plt.imshow(opt.mapWithObstacles)
	#plt.axis('tight')
	#plt.axis('equal')
	plt.show()
	r=10
	ang=np.arange(0,2*pi+0.01,0.01)
	x1=65
	y1=40
	x2=35
	y2=75
	plt.axis((opt.xmin,opt.xmax,opt.ymin,opt.ymax))
	
	
	plt.figure(1)
	ax1=plt.subplot(1,3,1)

	plt.imshow(np.flipud(opt.mapWithObstacles))
	#plt.axis('tight')
	#plt.axis('equal')
	plt.title('Candidate trajectories')

	ax3=plt.subplot(1,3,3)
	

	#check opt.htraj
	plt.title('Time-Average Statistics')
	plt.axis((opt.xmin,opt.xmax,opt.ymin,opt.ymax))
	
	agents.trajFigOPTIMAL=[0,0,0]
	agents.trajFig=[0,0,0]

	for iagent in range(0,opt.nagents):
		opt.iagent=iagent
		opt.figg=plt.figure(1)
		ax2.scatter(agents.xi[iagent,0],agents.xi[iagent,1],0,30)
		agents.trajFigOPTIMAL[iagent]=draw_path(traj(0*opt.z,opt,agents),opt.colors[iagent],2,5)

		#subplot

		agents.trajFig[iagent]=draw_path(traj(0*opt.z,opt,agents),opt.colors[iagent],2,5)
		ax1.scatter(agents.xi[iagent,0],agents.xi[iagent,1],0,30)

	Fig=draw_path(0*np.random.rand(3,3),'b',3,5)
	
	return opt,agents,Fig

	
		