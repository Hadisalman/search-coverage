#!/usr/bin/env python

import numpy as np
from numpy import matrix
import time
from initialization import initialization
from GenerateUtilityMap import GenerateUtilityMap
from GetFourierCoeff import GetFourierCoeff
from SMC_Update import SMC_Update
from Calculate_Ergodicity import Calculate_Ergodicity
from IPython import embed
import rospy
from visualization_msgs.msg import Marker
from visualization_msgs.msg import MarkerArray

class opt1(object):
    __slots__=('ergNkx','ergNky','ergmuk','ergHK','ergKX','ergKY','ergmu','ergs')
    __slots__=('Lx','Ly','nagents','Nsteps','dt','vlb','vub','wlb','wub','xmin','xmax','ymin','ymax')
class pos(object):
	__slots__=('x','y','theta')

#rviz visualization
#############################################################
topic='visualization_marker_array'
publisher=rospy.Publisher(topic,MarkerArray,queue_size=1000)
rospy.init_node('vis_marker',anonymous=True)
#rate = rospy.Rate(1000)
markerArray=MarkerArray()
marker1=Marker()
marker2=Marker()
marker3=Marker()
marker1.header.frame_id = "/my_frame"
marker1.id=0
marker1.type = marker1.SPHERE
marker1.action = marker1.ADD
marker1.scale.x = 10
marker1.scale.y = 10
marker1.scale.z = 10
marker1.color.a = 1.0
marker1.color.r = 0.0
marker1.color.g = 1.0
marker1.color.b = 0.0
marker2.header.frame_id = "/my_frame"
marker2.id=1
marker2.type = marker2.SPHERE
marker2.action = marker2.ADD
marker2.scale.x = 10
marker2.scale.y = 10
marker2.scale.z = 10
marker2.color.a = 1.0
marker2.color.r = 1.0
marker2.color.g = 0.0
marker2.color.b = 0.0
marker3.header.frame_id = "/my_frame"
marker3.id=2
marker3.type = marker3.SPHERE
marker3.action = marker3.ADD
marker3.scale.x = 10
marker3.scale.y = 10
marker3.scale.z = 10
marker3.color.a = 1.0
marker3.color.r = 1.0
marker3.color.g = 1.0
marker3.color.b = 0.0
#############################################################

if __name__ == "__main__":
	pose=pos
	opt=opt1
	
	opt,pose=initialization(opt,pose)
	addnoise=0

	X,Y,informationMap=GenerateUtilityMap(opt,addnoise)

	opt.ergmu=np.reshape(informationMap,X.shape)
	opt.ergmu=opt.ergmu/(opt.ergmu).sum()

	opt=GetFourierCoeff(opt,X,Y)

	dt=opt.dt
	Nsteps=opt.Nsteps
	Ck=matrix(np.zeros((opt.ergNkx,opt.ergNky)))
	Ergodicity_Metric=matrix(np.zeros((Nsteps,1)))
	
	#rviz visualization
	################################
	marker1.pose.orientation.w = 1.0
	marker1.pose.position.z = 20
	marker2.pose.orientation.w = 1.0
	marker2.pose.position.z = 20
	marker3.pose.orientation.w = 1.0
	marker3.pose.position.z = 20
	#################################
	
	t=time.time()
	for it in range(0,Nsteps):
		ti=(it+1)*dt
		pose,Ck=SMC_Update(pose,opt,ti,Ck)

		#rviz visualization
		########################
		marker1.pose.position.x = pose.x[0]
		marker1.pose.position.y = pose.y[0] 
		
		
		marker2.pose.position.x = pose.x[1]
		marker2.pose.position.y = pose.y[1] 
		
		
		marker3.pose.position.x = pose.x[2]
		marker3.pose.position.y = pose.y[2] 
		
		markerArray.markers.append(marker1)
		markerArray.markers.append(marker2)
		markerArray.markers.append(marker3)
		publisher.publish(markerArray)
		
		markerArray.markers.pop(0)
		markerArray.markers.pop(0)
		markerArray.markers.pop(0)
		########################
		ck=Ck/opt.nagents/ti
		Ergodicity_Metric=Calculate_Ergodicity(Ergodicity_Metric,ck,opt)
		#rate.sleep()
	
	print("It took:")
	print(time.time()-t)	





