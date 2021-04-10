import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate, optimize

#==============================================================================#
# UTILITY FUNCTIONS                                                            #
#==============================================================================#

def sig2a(pdf,lb,ub,eps,sig):
	""" For pdf 'pdf', this function returns the parameter 'a' so that the
		standard deviation of the random variable with pdf 'pdf' is 'sig'.

		pdf : a callable python function
			The pdf with arguments 'x', 'a', and 'eps'
		eps : float 
			The support parameter
		sig : float
			The desired standard deviation scaled by eps

	"""

	# standard deviation of the symmetric beta distribution with parameter 'a' 
	def std(a): 
		""" Computes the variance of the noise given parameter 'a'. 
		"""
		return np.sqrt(integrate.quad(lambda x: pdf(x,a,eps)*x**2,-eps,eps)[0])
	# solve for the parameter 'a' for which the noise has variance sig0**2
	print(std(lb)/eps)
	print(sig)
	print(std(ub)/eps)
	print("--------------")
	return optimize.bisect(lambda a: std(a)/eps-sig,lb,ub)
	#return optimize.newton(lambda a: std(a)/eps-sig,2.*eps)

sig2a = np.vectorize(sig2a) # vectorize 
		
def betapdf(x,a,eps):
	""" The beta distribution"""	
	return stats.beta.pdf(x,a,a,loc=-eps,scale=2*eps)

def normpdf(x,a,eps):
	"""the truncated normal distribution"""
	pdf_raw = stats.norm.pdf(x,loc=0,scale=a)
	cdf_dn  = stats.norm.cdf(-eps,loc=0.,scale=a)
	cdf_up  = stats.norm.cdf(+eps,loc=0.,scale=a)
	return pdf_raw/(cdf_up-cdf_dn)

#==============================================================================#
# THE ERROR DISTRIBUTION CLASS                                                 #
#==============================================================================#

class ErrorDistribution():

	def __init__(self,pdf,parasupp,distname,sig0=.95/np.sqrt(3),eps0=.2,m=3,n=20):

		# native parameters
		self.sig0 = sig0   # standard deviation of the noise
		self.eps0 = eps0   # the support of the noise is [-eps,+eps]

		# error function
		self.pdf = pdf # callable "error" function
		self.parasupp = parasupp # support of the distribution-specific parameter
		self.distname = distname

		# hard coded approximation parameter
		self.m = m # number of pdfs 
		self.n = n # number of standard deviations

		# sampling parameters: epsilon
		self.epsmax = .25 # analytical upper-bound
		self.epsmin = .05*self.epsmax

		# sampling parameters: sigma
		self.sigmax = .9/np.sqrt(3.) # analytical upper-bound
		self.sigmin = .5*self.sigmax

		# pdfs example standard deviations 
		self.sigex = [self.parasupp[0],(self.parasupp[0]+self.parasupp[1])/2,self.parasupp[1]]
		self.nvec = len(self.sigex)

	def cost_identity(self,eps,sig):
		""" computes the error of the identity message function"""
		# compute parameter for sigma
		print(self.parasupp)
		a = sig2a(self.pdf,self.parasupp[0],self.parasupp[1],eps,sig) 
		# middle
		def i_m(e): # integrand
			return self.pdf(e,a,eps)*e**2
		I_m = (1-2*eps)*integrate.quad(i_m,-eps,eps)[0] 
		# bottom-left
		def i_b(e,q):
			return ((q-(q+e+eps+0.)/2)**2)*self.pdf(e,a,eps)
		I_b = integrate.dblquad(i_b,0,2*eps,lambda q: -eps,lambda q: eps-q)[0]
		# top-right
		def i_t(e,q):
			return ((q-(q+e-eps+1.)/2)**2)*self.pdf(e,a,eps)
		I_t = integrate.dblquad(i_t,1-2*eps,1.,lambda q: 1.-eps-q,lambda q: eps)[0]
		return I_m+I_b+I_t

	# driver for discrete messages on the unit interval
	def cost_discrete(self,eps):
		""" computes the error of the identity message function"""
		return ((eps/(2*eps+1))**2)/3

#==============================================================================#

def plot_expdfs(ed,print2csv=True):

	nvec = 1000
	evec = np.linspace(-ed.eps0,ed.eps0,nvec)
	f = lambda e, a: ed.pdf(e,a,ed.eps0)
	for j in range(1,len(ed.sigex)):
		plt.plot(evec,f(evec,ed.sigex[j]))
	plt.plot(0.,0.)
	plt.grid()
	#plt.ylim([0,])
	plt.xlabel('error')
	plt.title(ed.distname + ' PDFs')
	plt.show()

	# print example PDFs
	if print2csv:

		# example PDFs
		evec = np.linspace(-ed.eps0,ed.eps0,ed.nvec)
		f = lambda e, a: ed.pdf(e,a,ed.eps0)
		with open(ed.distname + '_pdfs.csv','w') as csvfile:
			pen = csv.writer(csvfile)
			pen.writerow(['e','PDF1','PDF2','PDF3'])
			for j in range(0,ed.nvec):
				row = np.hstack([evec[j],f(evec[j],ed.sigex[0]),f(evec[j],ed.sigex[1]),f(evec[j],ed.sigex[2])])
				pen.writerow(row)


def plot_errors(ed,dim,print2csv=True):
	""" This function performs three functions. First, it computes the error
		associated with the parameters varied along 'dim'. If 'dim' = 'sig',
		for example, then 'eps' is held fixed, and the error of the discrete
		and identity message functions is computed. Second, it plots the
		result. Third (and optionally), it writes the results to a csv file.

		ed : ErrorDistribution
			The error distribution to plot
		dim : str
			The dimension (parameter) to vary (either 'eps' or 'sig')
		print2csv : logical
			If true, then print the results to .csv file for external use

	"""

	# vectorize the cost functions
	cost_identity = np.vectorize(ed.cost_identity)
	cost_discrete = np.vectorize(ed.cost_discrete)
		
	def int_plot(epsvec,sigvec,dim,xlabel,title):
		"""This is the internal plotting function"""

		# compute costs
		CD = cost_discrete(epsvec)
		CI = cost_identity(epsvec,sigvec)

		# plot
		if dim == 'eps':
			plt.plot(epsvec,CD)
			plt.plot(epsvec,CI)
		else:		
			plt.plot(sigvec,CD)
			plt.plot(sigvec,CI)
		plt.legend(['discrete','identity'])
		plt.xlabel(xlabel)
		plt.ylabel('expected loss')
		plt.title(title)			
		plt.grid()
		plt.show()
		return CD, CI

	# hold sigma fixed, vary epsilon
	if dim == 'eps': 

		epsvec = np.linspace(ed.epsmin,ed.epsmax,ed.n)
		sigvec = ed.sig0*np.ones(ed.n)

		title = 'support of noise: ' + ed.distname
		CD, CI = int_plot(epsvec,sigvec,dim,'$\\bar{\\epsilon}$',title)

	# hold epsilon fixed, vary sigma
	elif dim == 'sig': 

		epsvec = ed.eps0*np.ones(ed.n)
		sigvec = np.linspace(ed.sigmin,ed.sigmax,ed.n)

		xlabel = '$\\sigma/\\bar{\\epsilon}$'
		title = 'scaled standard deviation of noise: ' + ed.distname
		CD, CI = int_plot(epsvec,sigvec,dim,xlabel,title)

	else:

		print("warning: unrecognized dimension")

	# print example PDFs
	if print2csv:

	    # costs of discrete and identity message functions
		filename = ed.distname + '_' + dim + '_discrete_identity.csv' 
		with open(filename,'w') as csvfile:
			pen = csv.writer(csvfile)
			pen.writerow([dim,'cost_discrete','cost_identity'])
			for j in range(0,ed.n):
				if dim == 'eps':
					pen.writerow(np.hstack([epsvec[j],CD[j],CI[j]]))
				else:
					pen.writerow(np.hstack([sigvec[j],CD[j],CI[j]]))