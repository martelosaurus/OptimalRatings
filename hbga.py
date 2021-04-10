"""This is a home brewed genetic algorithm (hbga) module"""
import numpy as np

# parameters
n = 100 # number of test points

#==============================================================================#
#   CLASSES                                                                    #
#==============================================================================#
class Model:
"""a model contains specifications for the genetic algorithm"""

	def __init__(self,val):

		# native
		self.n = 1000 # number of points in the partition (including 0. and 1.)
		self.maxiter = 1000 # maximum number of iterations
		self.popsize = 5000 # size of the initial population
		self.numchld = 2    # number of children
		self.killpct = .5   # percent of the population to kill

		# procreation parameters
		self.w = .1 # weight placed on mutation

#==============================================================================#
#   UTILITY FUNCTIONS                                                          #
#==============================================================================#

def randpart(n):
"""returns a random partition of [0,1] of size n"""
	x = np.random.rand(n)     # n random points in [0,1]
	x = np.cumsum(x)          # create an increasing sequence using cumsum
	x = (x-x[0])/(x[-1]-x[0]) # scale so that x[0] = 0 and x[-1] = 1
	return x

def chldpart(x,y,n):
"""returns a child partitions"""
	genes = np.vstack([x,y,randpart(n)])
	return np.mean(gene,0,[1-np.w,1-np.w,np.w])

#==============================================================================#
#   LOSS FUNCTION                                                              #
#==============================================================================#

def actfunc(m_tilde):
"""returns the receiver's optimal action having received m_tilde"""
	q_dn = 1. 
	q_up = 2.

def lossfun(m):


#==============================================================================#
#   EVOLUTION                                                                  #
#==============================================================================#

def evol(,model):
"""evolves a random population of solution with the specifications in 'model'"""

	population = randpart(model.popsize)

	while True:
		# step 1: kill weak solutions


		# step 2: procreate
		population.procreate()



	return minpart





