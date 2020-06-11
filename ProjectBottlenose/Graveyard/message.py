"""This is a home brewed genetic algorithm (hbga) module"""
"""   Parents have enough children to keep the population size constant"""

import numpy as np
from scipy import linalg, integrate, optimize, stats
import matplotlib.pyplot as plt

# default parameter values
eps0 = .1
sig0 = .05

# variance of the symmetric beta distribution with parameter 'a' 
def BetaVar(a): 
	return stats.beta.expect(lambda x: x**2,(a,a),-eps0,2*eps0) 

# solve for the parameter 'a' for which the noise has variance sig0**2
a0 = optimize.newton(lambda a: BetaVar(a)-sig0**2,x0=2.)

# default importance function: uniform
def iFun0(x):
	return stats.beta.pdf(x,1.,1.) 

# default error function: beta
def eFun0(x):
	return stats.beta.pdf(x,a0,a0,-eps0,2*eps0) 

class Message():
	""" A message function class
		
		Notes
		-----
		To save on computation, eps is rounded so that the distance between 
		knots evenly divides eps. Therefore, if m is an element of the set of 
		knots M, then m+eps and m-eps are also elements of M.  

	"""	


	def __init__(self,sig=sig0,eps=eps0,M_dn=0.,M_up=1.,iFun=iFun0,eFun=eFun0):

		# native parameters
		self.sig = sig   # standard deviation of the noise
		self.eps = eps   # the support of the noise is [-eps,+eps]
		self.M_dn = M_dn # lower bound of the message space
		self.M_up = M_up # upper bound of the message space

		# error and importance functions
		self.eFun = eFun # callable "error" function
		self.iFun = iFun # callable "importance" function 

		# hard coded approximation parameters
		# NOTE: the parition includes 0. but not 1.
		self.n = 1000                    # parition size (with 0. but not 1.)
		self.h = 1./self.n                      # partition mesh
		self.M = np.arange(0.,1.+self.h,self.h) # the parition (with 0. and 1.) 
		self.L = self.M

		# run parameters
		self.max = 1000  # maximum number of iterations
		self.tol = 1.e-10 # error tolerance
		self.eqtol = 1.e-10

		# round eps
		self.eps = (np.floor(self.eps*self.n)/self.n)
		self.EPS = np.arange(-self.eps,self.eps+self.h,self.h)


	def Plot(self):
		"""Plots the message function"""
		plt.plot(self.M,self.M,'-b')
		plt.plot(self.L,self.M,'--m')
		plt.axis([0.,1.,self.M_dn,self.M_up])
		plt.xlabel('state $(q)$')
		plt.ylabel('message $(m)$')
		plt.grid()
		plt.show()

	def Optimize(self):

		# initialize the error
		error = 2*self.tol

		# guess the solution
		L = self.M

		# contract
		it = 0
		while (error>self.tol) and (it<self.max):

			print(it)

			L_ref = L

			it += 1

			def action():		
				neps = np.int(2*self.eps*self.n)
				q_dn = L[np.clip(np.r_[0:self.n+1]-neps,0,None)] 
				q_up = L[np.clip(np.r_[0:self.n+1]+neps,None,self.n)]
				A_up = (q_up+L)/2
				A_dn = (q_dn+L)/2

				return A_up, A_dn, 

			# solve for the next iterate by solving the quadratic
			A_up, A_dn = action()
			L = .5*(A_up+A_dn)
			L[0] = 0.
			L[-1] = 1.

				# compute the error
			if it%100==0:
				plt.plot(L)
				plt.show()

			error = np.linalg.norm(L-L_ref)

			#break
		self.L = L