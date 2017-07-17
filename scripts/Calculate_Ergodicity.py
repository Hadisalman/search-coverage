import numpy as np
from numpy import matrix


def Calculate_Ergodicity(Ergodicity_Metric,ck,opt):

	LK=opt.ergLK
	muk=opt.ergmuk

	Ergodicity_Metric=(matrix(np.multiply(LK,np.square(ck-muk)))).sum()

	return Ergodicity_Metric


