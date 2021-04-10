"""This is a home brewed genetic algorithm (hbga) module"""
"""   Parents have enough children to keep the population size constant"""

import numpy as np
from scipy import integrate, optimize, stats
import genetic
import matplotlib.pyplot as plt

#==============================================================================#
#   MESSAGE CLASS                                                              #
#==============================================================================#

# default parameter values
eps0 = .1
sig0 = .05

# variance of the symmetric beta distribution with parameter 'a' 
var = lambda a: stats.beta.expect(lambda x: x**2,(a,a),-eps0,2*eps0) 

# solve for the parameter 'a' for which the noise has variance sig0**2
a0 = optimize.newton(lambda a: var(a)-sig0**2,x0=2.)
a0 = 1.

# default importance function: uniform
def iFun0(x):
	return stats.beta.pdf(x,1.,1.,-eps0,2*eps0) 

# default error function: beta
def eFun0(x):
	return stats.beta.pdf(x,a0,a0,-eps0,2*eps0) 

def Moment(Fun,k,A,B):
	""" Integrates x^k against iFun(x) on [A,B] using the trapezoid rule 
	
		Parameters
		----------
		Fun : callable Python function
			Either self.eFun or self.iFun.
		k : float
			Moment.
		A : ndarray of dim (?,partsize)
			Lower integration bounds.
		B : ndarray of dim (?,partsize)
			Upper integration bounds.
	
		Notes
		-----
		Because [A,B] is typically small relative to [0,1], this function uses a
		simple trapezoid rule to compute the integral.
	"""
	return .5*(B-A)*((B**k)*Fun(B)-(A**k)*Fun(A))


class Message():
	"""A message function class"""

	def __init__(self,sig=sig0,eps=eps0,M_dn=0.,M_up=1.,iFun=iFun0,eFun=eFun0):

		# native parameters
		self.sig = sig   # standard deviation of the noise
		self.eps = eps   # the support of the noise is [-eps,+eps]
		self.M_dn = M_dn # lower bound of the message space
		self.M_up = M_up # upper bound of the message space

		# error and importance functions
		self.eFun = eFun # callable "error" function
		self.iFun = iFun # callable "importance" function 

		# hard coded genetic algorithm parameters
		self.partsize = 100 # size of the partition (including 0. and 1.)
		self.maxiter  = 1000  # maximum number of iterations
		self.popsize  = 10000 # size of the initial population
		self.numchld  = 2    # number of children
		self.killpct  = .5   # percent of the population to kill

		# pre-compute the loss functions
		self.q  = np.linspace(0.,1.,self.partsize) # partition in q

		# pre-allocate space for iMoments
		self.Iq = [np.zeros((self.popsize,self.popsize)) for j in range(0,3)]

		# procreation parameters
		self.w = .75 # weight placed on mutation

	def Optimize(self):
		"""Wrapper for genetic.evolution"""	
		self.m, self.loss = genetic.evolution(self)
		return self.m, self.loss

	def Plot(self):
		"""Plots the message function"""
		plt.plot(self.q,self.m,'-b')
		plt.plot(self.q,self.q,'--m')
		plt.axis([0.,1.,self.M_dn,self.M_up])
		plt.xlabel('state $(q)$')
		plt.ylabel('message $(m)$')
		plt.grid()
		plt.show()

	def eMoment(self,k,a,b):
		""" Wrapper for Moment with Fun = self.eFun"""
		return Moment(self.eFun,k,a,b)

	def iMoment(self,k,a,b):
		""" Wrapper for Moment with Fun = self.iFun"""
		return Moment(self.iFun,k,a,b)


