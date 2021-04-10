"""Classes for Clarifying by Discretizing"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy import integrate

#==============================================================================#
#   PARAMETERS                                                                 #
#==============================================================================#
def Optimize(Model):
	"""This is a wrapper for HDBF, a "High Dimensional Brute Force" optimizer"""

#==============================================================================#
# CLASS: MESSAGE FUNCTION                                                      #
#==============================================================================#
class MessageFunction():

	def __init__(self,sig,eps,pdf,M_dn,M_up,iFun):

		# native parameters
		self.sig = sig   # standard deviation of the noise
		self.eps = eps   # the support of the noise is [-eps,+eps]
		self.pdf = pdf   # either 'beta' or 'normal'
		self.M_dn = M_dn # lower bound of the message space
		self.M_up = M_up # upper bound of the message space
		self.iFun = iFun # callable "importance" function 

		# optimization parameters
		N = 1000 # number of steps in the approximation

		# message functions
		self.Q = hbga.evolution()
		self.M = np.zeros()

		# message properties
		self.cost = 0.

	def plot(self):
		"""Plots the message function"""
		plt.plot(Q,M,'-b')
		plt.plot(M,M,'--m')
		plt.axis([M_dn-eps,M_up+eps,-.1,1.1])
		plt.xlabel('state $(q)$')
		plt.ylabel('message $(m)$')
		plt.grid()
		plt.show()

#==============================================================================#
# CHILD CLASS: DISCRETE MESSAGE FUNCTION                                       #
#==============================================================================#
class DiscreteMessageFunction(MessageFunction):

	def __init__(self,sig,eps,pdf,M_dn,M_up):
		super().__init__(self,sig,eps,pdf,M_dn,M_up)


#==============================================================================#
# CHILD CLASS: IDENTITY MESSAGE FUNCTION                                       #
#==============================================================================#
class IdentityMessageFunction(MessageFunction):

	def __init__(self,sig,eps,pdf,M_dn,M_up):
		super().__init__(self,sig,eps,pdf,M_dn,M_up)