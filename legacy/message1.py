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
a0 = 1.
def eFun0(x):
	#return stats.norm.pdf(x,0.,eps0/2)
	return stats.beta.pdf(x,a0,a0,-eps0,2*eps0) 
		
class Message():

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
		self.tol = 1.e-4 # error tolerance

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

			# moments of the importance function
			def moment(self,k,a,b):
				moms = np.zeros(len(a))
				for j in range(0,len(a)):
					#if a[j]>=b[j]:
						#print("warning: danger")
					I = lambda q: (q**k)*self.iFun(q)
					moms[j] = integrate.quad(I,a[j],b[j])[0]

				return moms

			# return the optimal action
			def action():		
				# largest and smallest possible states of the world
				neps = np.int(2*self.eps*self.n) # how many steps in
				q_dn = L[np.clip(np.r_[0:self.n+1]-neps,0,None)] 
				q_up = L[np.clip(np.r_[0:self.n+1]+neps,None,self.n)]
				#print(np.vstack([q_dn,L,q_up]).T)
				I1_dn = moment(self,1.,q_dn,L)
				I0_dn = moment(self,0.,q_dn,L)
				I1_up = moment(self,1.,L,q_up)
				I0_up = moment(self,0.,L,q_up)
				A_up = I1_up/I0_up
				A_dn = I1_dn/I0_dn
				A_dn[0] = 0.
				A_up[-1] = 1.
				
				#A_up_check = (q_up+L)/2
				#A_dn_check = (q_dn+L)/2
				#A_dn[0] = A_dn_check[0]
				#A_up[-1] = A_up_check[-1]
				#print(linalg.norm(A_up[:-1]-A_up_check[:-1]))
				#print([A_up[-1],A_up_check[-1]])
				#print(linalg.norm(A_dn[1:]-A_dn_check[1:]))
				#print([A_dn[0],A_dn_check[0]])
				# other set
				neps = np.int(self.eps*self.n) # how many steps in
				q_dn = L[np.clip(np.r_[0:self.n+1]-neps,0,None)] 
				q_up = L[np.clip(np.r_[0:self.n+1]+neps,None,self.n)]
				#print(np.vstack([q_dn,L,q_up]).T)
				I1_md = moment(self,1.,q_dn,q_up)
				I0_md = moment(self,0.,q_dn,q_up)
				A_md = I1_md/I0_md
				A_md = I1_md
				A_md[0] = 0.
				A_md[-1] = 1.

				return A_md, A_up, A_dn, 

			# remainder integral
			def remain(k,A_md):
				F = self.iFun(np.add.outer(self.M,self.EPS))
				A = linalg.hankel(A_md,np.flipud(A_md[-len(self.EPS):])) 
				#print(A)
				return np.sum(.5*(F[:,1:]-F[:,:-1])*(A[:,1:]**k+A[:,:-1]**k),1) 

			# solve for the next iterate by solving the quadratic
			A, A_up, A_dn = action()
			T1 = A_up*self.eFun(self.eps)-A_dn*self.eFun(-self.eps)
			T2 = remain(0.,A)*(L**2)-2*remain(1.,A)*L+remain(2.,A)
			L = -.5*(T2/T1-(A_up+A_dn))
			L = (A_up+A_dn)/2
			# compute the error
			error = np.linalg.norm(L-L_ref)
			#print(error)

		# set the (inverse message function)
		self.L = L