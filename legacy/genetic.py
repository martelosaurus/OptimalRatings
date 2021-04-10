"""This is a home brewed genetic algorithm module"""

import numpy as np
from scipy import integrate, optimize

def randpart(m,n):
	"""	Help on function actfun in module genetic
		
		Parameters
		----------
		m : int 
			The number of partitions to return. 
		n : int
			The size of each partition (including 0. and 1.)

		Returns
		-------
		x : ndarray
			An ndarry of shape (m,n), where each row is a partition.

		Notes
		-----
		A partitions is constructed in three steps: n points are
			1. drawn according to a uniform distribution on [0,1]
			2. mapped to an increasing sequence of numbers using cumsum
			3. and scalled so that x[0] = 0 and x[n-1] = 1
	"""

	x = np.random.rand(m,n)      
	x = np.cumsum(x,1) 
	for j in range(0,n):
		x[:,j] = (x[:,j]-x[:,0])/(x[:,n-1]-x[:,0]) 
	return x

def procreate(fit,msg):
	"""	Help on function actfun in module genetic
		Returns child partitions from a set of parent partitions.

			Parameters
	      	---------- 
	        fit : np.array of dim (popsize/2,partsize)
				An ndarray of fit parent partitions.

			Returns
			-------
			m_children : ndarray of dim (popsize,partsize)
				An ndarray of children partitions.

			Notes
			-----
			The function chooses random pairs of surviving parents to produce 
			four children. 
	"""

	# mothers
	mothers = fit[0:msg.popsize-1:2,:] # even indices
	mothers = np.tile(mothers,(4,1))

	# fathers
	fathers = fit[1:msg.popsize:2,:] # odd indices
	fathers = np.tile(fathers,(4,1))

	# zombies
	zombies = randpart(msg.popsize,msg.partsize) 


	np.tile(np.average(fit,0),(msg.popsize,1))
	return .5*msg.w*mothers+.5*msg.w*fathers+(1-msg.w)*zombies

def actfunc(m,m_tilde,msg):
	"""	Help on function actfun in module genetic
	 	Returns the receiver's optimal action having received m_tilde.
			
			Parameters
	      	---------- 
	        m : np.array
				A partition.
			m_tilde : float
				A received message on which the receiver acts.

			Returns
			-------
			action : float
				The receiver's optimal action having received m_tilde.

			Notes
			-----
			l(m_tilde) is an interval in [0,1]. Take the inf for q_dn and the 
			sup for q_up.
	"""

    # largest and smallest possible states of the world
	q_tilde = q[(m-eps<=m_tilde) & (m_tilde <= m+eps)]
	q_dn = np.min(q_tilde)
	q_up = np.max(q_tilde)

	# moments of the error distribution
	I = lambda k: integrate.fixed_quad(lambda q: (q**k)*model.iFun(q),q_dn,q_up) 

	# return the optimal action
	return I(1.)/I(0.)

def fitness(M,msg):
	"""	Help on function actfun in module genetic
		Computes the loss for a population of partitions.

			Parameters
	      	---------- 
	        M : ndarray of dim (popsize,partsize)
				An ndarray of partitions.
			msg : a Message object
				A message object containing specifications.

			Retruns
			-------
			loss : float
				The expected loss incurred by using the message function m. 	

			Notes
			-----
			The function computes the loss in each of the partsize-1 state bins,
			[j/n,(j+1)/n]. On each such bin, the message function is constant 
			and the loss can be computed with one-dimensional definite 
			integrals. Because the action is dicontinuous at m_tilde = 2*eps and
			m_tilde = 1-2*eps, the loss is computed separately in each of five
			regions in q-e space. 
	"""

	#M = M[:,:-1] # use left step-function
	M = M[:,1:] # use right step-function

	def LossTerm(qk,ek):
		""" Loss Terms
			
			Parameters
			----------

			qk : float
				q-moment.
			ek : e-moment

			Notes
			-----
			The state 'q' lies on the x-axis and the noise 'e' on the y-axis. 
			"left", "middle", and "right" describe the x-axis, while "up" and 
			"down" describe the y-axis. 

		"""

		# pre-allocate
		e_md = np.zeros((msg.popsize,msg.partsize-1))		
		
		# middle
		e_md[M<=2*msg.eps] = 2*msg.eps-M[M<=2*msg.eps]
		e_md[1-2*msg.eps<=M] = 1-2*msg.eps-M[1-2*msg.eps<=M]

		# don't hate...integrate
		e_dn_int = msg.eMoment(ek,-msg.eps,e_md)
		e_up_int = msg.eMoment(ek,e_md,+msg.eps)

		return msg.Iq[np.int(qk)]*(e_up_int+e_dn_int)

	# compute the loss on each of the five continuous regions in q-e space
	loss = LossTerm(2.,0.)-2.*LossTerm(1.,1.)+LossTerm(0.,2.)

	# sum across the losses over the state bins
	loss = np.sum(loss,1) 

	# fitness = -loss
	return -loss

def evolution(msg):
	""" Help on function evoluation in module genetic
	    Evolves a random population of parition with the specifications in model 
	    the population evolves for 'maxiter' periods and then the fittest 
	    parition is returned.

		Parameters
		---------- 
	    msg : A Message object
			A Message object containing model specifications.

		Return
		------
		most_fit_part : ndarray
			The partition associated with the lowest loss.
		most_fit_part : float
			The loss associated with the most fit parition.
	""" 
	# precompute loss integrals on each of the partsize-1 bins
	for k in range(0,3):
		temp = msg.iMoment(1.*k,msg.q[:-1],msg.q[1:])
		msg.Iq[k] = np.tile(temp,(msg.popsize,1))

	# store population in a matrix (each row is a parition)
	population = randpart(msg.popsize,msg.partsize)

	# iterate
	for iter in range(0,msg.maxiter):

		print(iter)

		# step 1: find the fit parents
		fit = fitness(population,msg)
		I = np.argsort(-fit)
		I = I[:msg.popsize/2] # kill the bottom half

		# step 2: fit parents give birth to the next generation
		if iter<msg.maxiter-1:
			population = procreate(population[I,:],msg)

	# return the most fit partition and its associated loss
	mostfitpart = population[I[0],:] # partition
	mostfitloss = -fit[I[0]]         # loss = -fitness
	return mostfitpart, mostfitloss 