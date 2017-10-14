# ///-----------------------------------------------------------------
# ///   Description:    <STOEC>
# ///   Author:         <Hadi Salman>                    
# ///   Date:           <Oct 13th 2017>
# ///   Revision History: ---
# ///-----------------------------------------------------------------
import numpy as np
from math import exp,pi,sqrt,log
from collections import namedtuple
import scipy
from IPython import embed
from rmvnrnd import rmvnrnd

def entropy(mun,S):
	k=len(mu)
	f=k/2*(1+log(2*pi))+np.log(np.linalg,det(S))/2
	return f

def minb3(b,S):
	N=S.Js.size
	Jmin=np.min(S.Js)
	Jmax=np.max(S.Js)

	ws=np.exp(-b*S.Js/Jmax)

	delta=0.05
	g=sqrt(log(1/delta)/(2*N))

	f=log(np.mean(ws)-exp(-b*Jmin/Jmax)*g)+b*(np.mean(S.Js)/Jmax - g)
	return -f


def cem(fun,x0,opts,varargin,agents,erg):

	S=namedtuple('S',['Js','Jh','ps0','ps','mu','xs','v','C'])
	
	d=x0.shape[0]

	if hasattr(opts,'N')==False:
		opts.N=d*20
	if hasattr(opts,'rho')==False:
		opts.rho=0.1
	if hasattr(opts,'C')==False:
		opts.C=np.identity(d)
	if hasattr(opts,'iter')==False:
		fun=20
	if hasattr(opts,'v')==False:
		opts.v=1
	if hasattr(opts,'sigma'):
		N=2*d+1
	if hasattr(opts,'tilt')==False:
	 	opts.tilt=0
	if hasattr(opts,'lb')==False:
		opts.lb=[]
	if hasattr(opts,'ub')==False:	
	 	opts.ub=[]
	 	
	N=opts.N
	nf=round(opts.rho*N)
	C=opts.C
	v=opts.v

	cs=np.zeros((N,1))
	xs=np.zeros((d,N))
	x=x0
	c=float("inf")

	mu=np.array(x0).T[-1]

	a=0.001
	k=0
	b=2
	l=a*a*(d+k)-d

	Ws= [l/(d+l)] + [ 1/(2*(d+l)) for i in range(2*d)]
	Wc= [l/(d+l) + (1-a*a+b)] + [1/(2*(d+l)) for i in range(2*d)]

	for j in range(1,opts.iter+1):
		
		# if hasattr(opts,'sigma'):
		# 	A=sqrt(d+l)*np.linalg.cholesky(C).T
		# 	xs= [] 
			
		# 	xm=np.zeros((d,1))
		# 	for i in range(0,xs.shape[1]):
		# 		fi=fun(xs[:,i],varargin)
		# 		if len(fi)>1:
		# 			cs[i,0]=np.multiply(fi,fi).sum()
		# 		else:
		# 			cs[i,0]=fi
		# 		cs[i,0]=exp(-cs[i,0])
		# 		xm=xm +Ws[i]*cs[i,0]*xs[:,i]

		# 	Pm=np.zeros((d,d))
		# 	for i in range(0,xs.shape[1]):
		# 		dx=np.matrix(xs[:,i]-xm)
		# 		Pm=Pm+Wc[i]*cs[i,0]*dx*dx.T

		# 	csn=cs.sum()
		# 	mu=mu/csn
		# 	C=Pm/csn
		# 	x=mu
		# 	c=cs[0,0]
		# else:
		if  opts.lb.T.tolist()[0]:
			n=mu.size
			A=np.concatenate((-np.identity(n),np.identity(n)))
			B=np.concatenate((-opts.lb,opts.ub))
			# embed()
			xs=rmvnrnd(mu,C,N,A,B,100,0)
			# print(xs[20])
		else:
			xs=np.random.multivariate_normal(mu,C,N)
			xs=xs.T
		xs=xs.T
		# embed()
		for i in range(0,N):
			
			fi=fun(xs[:,i],varargin,agents,erg)
			
			if len([fi])>1:
				cs[i]=fi*fi/2
			else:
				cs[i]=fi

		if not opts.tilt:

			iss=np.argsort(cs.flatten())
			cs= np.sort(cs.flatten())
			
			xes=xs[:,iss.flatten().tolist()[0:int(nf)]]

			xes =np.reshape(xes,xes.shape[0:2],'F')

			mu=(1-v)*mu[-1]+v*np.mean(xes.T,axis=0)
			# embed()
			C=np.multiply((1-v),C)+v*np.cov(xes)
			# embed()
			if cs[0]<c:
				x=xes[:,0]
				c=cs[0]

		else:
			if j==1:
				mvn1 = multivariate_normal(mu.T,C)
				S.ps0 = mvn1.pdf(xs.T)

			imin=np.argmin(cs)
			cmin=np.min(cs)
			b=1/np.mean(cs)	

			# if False:
			# 	b=b*(entropy(mu,C))

			# 	mvn2 = multivariate_normal(mu.T,C)
			# 	S.ps = mvn2.pdf(xs.T)
			# 	S.Jh=np.mean(cs)
			# 	S.Js=cs
			# 	bmin=0
			# 	bmax=1
			# 	S.xs=xs
			# 	S.v=v
			# 	S.mu=mu
			# 	S.C=C
			# 	bs=np.arange(bmin,bmax+0.001,0.001)
			# 	gs=np.zeros(bs.shape)
			# 	for l in range(0,len(bs)):
			# 		gs(l)=minb3(bs[l],S)

			# 	#plotting

			# 	bi=np.argmin(gs)
			# 	gm=np.min(gs)
			# 	b=bs[bi]
			# 	b,FVAL,EXITFLAG,OUTPUT=scipy.optimize.fminbound(minb3,bmin,bmax)

			ws=np.exp(-b*cs)
			ws=ws/ws.sum()

			mu=(1-v)*mu+v*(xs*ws)
			C=(1-v)*C	+v*weightedcov(xs.T,ws)
			if cmin<c:
				x=xs[:,imin]
				
				c=cmin

	return x,c,mu,C				



# def weightedcov(Y,w):

# 	ctrl=isinstance(w,ndarray) & np.imag(w)==0 & not np.nan(w) & not np.inf(w)






					