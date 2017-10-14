# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------
from math import pow
import numpy as np
def Calculate_Ergodicity(ck,muk,Nagents):


	Nkx=muk.shape[0]
	Nky=muk.shape[1]
	Ergodicity_Metric=0
	for iagent in range(0,Nagents):
		for kx in range(0,Nkx):
			for ky in range(0,Nky):
				lambda_k=1.0/pow(1.0+kx*kx+ky*ky,1.5)
				Ergodicity_Metric=Ergodicity_Metric+lambda_k*pow(abs(ck[kx,ky]-muk[kx,ky]),2)


	return Ergodicity_Metric				