# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def draw_path(xs,c,lw,ms):

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	G = x.plot(xs[0,:], xs[1,:], xs[2,:])
	plt.setp(G,linewidth=lw,markersize=ms)

	return G		