import numpy as np
from math import pi
import scipy.special
import scipy.optimize
from IPython import embed


def asymcdfinv(y,erf1,de,a):

	f=erf1+de*y
	l=np.log(np.multiply(-f,2+f))
	b=2/(pi*a)+l/2
	x=np.sqrt(-b+np.sqrt(np.square(b)-l/a))

	return x

def meantrunc(lower,upper,s):

	if s ==float("inf"):
		m=(upper+lower)/2
		
	else:
		a=(lower/np.sqrt(2))/s
		b=(upper/np.sqrt(2))/s
		
		corr=np.sqrt(2/pi)*np.divide((-np.exp(-np.square(b))),scipy.special.erf(b)-scipy.special.erf(a))
		m=np.multiply(s,corr)

	return m

def vartrunc(lower,upper,s):
	
	if s == float("inf"):
		v=pow((upper-lower),2/12)
	else:
		a=(lower/np.sqrt(2))/s
		b=(upper/np.sqrt(2))/s
		if a ==float("inf"):
			ea=0
		else:	
			ea=np.multiply(a,np.exp(-np.square(a)))
				
		if b == float("inf"):
			eb=0
		else:
			eb=np.multiply(b,np.exp(-np.square(b)))

		corr=1-np.divide((2/np.sqrt(pi))*(eb-ea),scipy.special.erf(b)-scipy.special.erf(a))
		v= np.multiply(np.square(s),corr)	
			
	return v

def stdtrunc(lower,upper,s):
	
	temp1=vartrunc(lower,upper,s)
	temp2=meantrunc(lower,upper,s)
	arg=temp1-np.square(temp2)
	
	stdt=np.sqrt(np.abs(arg))

	return stdt

def scz(sc,*data):
	targetsigma2,lower,upper=data
	# embed()
	# print("hi")
	
	res=vartrunc(lower,upper,sc)-targetsigma2- np.square(meantrunc(lower,upper,sc))
	
	return res



def TruncatedGaussian(sigma,range,varargin):
	
	PREVSIGMA=sigma
	PREVRANGE=range
	PREVSIGMAC=0
	shapeflag=(sigma<0)
	range=np.float32(range).tolist()

	if isinstance(range,list)==False:
		range=abs(range)
		range=[-range,range]
	else:
		range=np.sort(range).tolist()

	sigma=np.float32(np.abs(sigma))
	
	
	n=varargin
	if sigma<0:
		sigmac=sigma
	else:
		
		if (np.square(np.diff(range)).tolist())[0]<12*pow(sigma,2):
			# embed()
			print("TruncatedGaussian:RangeSigmaIncompatible")
			sigmac=float("inf")
		
		elif [sigma,range]==[PREVSIGMA,PREVRANGE]:
			sigmac=PREVSIGMAC
			
		else:
			# embed()
			data=(pow(sigma,2),range[0],range[1])
			sigmac=scipy.optimize.fsolve(scz,sigma,args=data)
			sigmac=np.abs(sigmac)[0]
			# embed()
			# if flag<0:
				# print("error:TruncatedGaussian:fzerofailled")

			PREVSIGMA=sigmac
			PREVRANGE=range
			PREVSIGMAC=sigmac

	meaneffective=meantrunc(range[0],range[1],sigmac)
	sigmaeffective=stdtrunc(range[0],range[1],sigmac)

	if sigmac ==float("inf"):
		#if any
		cdfinv=lambda y: range[0]*y

	else:

		c=np.sqrt(2)*sigmac
		rn=np.divide(range,c)
		asymthreshold=4
		#if any
		if np.prod(np.sign(rn))>0 and np.abs(rn)>= asymthreshold:

			c=c*np.sign(rn[0])
			rn=np.abs(rn)
			left=np.min(rn)
			right=np.max(rn)

			a=0.147

			x2=left*right
			ax2=a*x2
			e1=(4/pi+ax2)/(1+ax2)
			e1=np.exp(-x2*e1)

			x2=right*right
			ax2=a*ax2
			e2=(4/pi+ax2)/(1+ax2)
			e2=np.exp(-x2*e1)

			de=-0.5*(e2-e1)-0.125*(e2-e1)*(e1+e2)

			erf1=(-0.5*e1-0.125*pow(e1,2))

			cdfinv=lambda y:np.multiply(asymcdfinv(y,erf1,de,a),c)
		else:
			
			e=scipy.special.erf(np.divide(range,c))	
			
			cdfinv=lambda y: np.multiply(scipy.special.erfinv(e[0]+np.multiply(np.diff(e),y)),c)
	# embed()		
	X=cdfinv(np.random.rand(n,n))	

	X=np.maximum(np.minimum(X,range[1]),range[0])	
	# embed()
	
	return X