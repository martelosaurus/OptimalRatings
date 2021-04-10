import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate
from joblib import Parallel, delayed
import multiprocessing

# driver for identity messages on the unit interval
def int_identity(eps,bsp):

	# action function
	def A(m_tilde): 
		q_dn = max(m_tilde-eps,0)
		q_up = min(m_tilde+eps,1)
		f = lambda e: m_tilde-e
		if q_up>q_dn:
			return stats.beta.expect(f,bsp,q_dn,q_up-q_dn,q_dn,q_up,True)
		else:
			return q_up
	A = np.vectorize(A)

	# loss function
	L = lambda e, q: (q-A(q+e))**2

	# boundary functions
	b1 = lambda q: eps-q
	b2 = lambda q: 1-eps-q

    # return the cost
    #   note: break the integral into five regions to resolve the discontinuity
    #         in the action function
	return (integrate.dblquad(L,0,2*eps,-eps,b1)[0]
		+integrate.dblquad(L,0,2*eps,b1,eps)[0]
		+integrate.dblquad(L,2*eps,1-2*eps,-eps,eps)[0]
		+integrate.dblquad(L,1-2*eps,1,-eps,b2)[0]
		+integrate.dblquad(L,1-2*eps,1,b2,eps)[0])

# driver for discrete messages on the unit interval
def int_discrete(eps,bsp):
	# The per-bin loss is just the error variance
	return stats(lambda x: x**2,bsp,-eps,2*eps,-eps,eps,True)

# driver
def processInput(j,N,eps_vec,bsp_vec):
	cost = np.zeros(N)
	for i in range(0,N):
		cost[i] = int_identity(eps_vec[i],(bsp_vec[j],bsp_vec[j])) 
	return cost

if __name__ == '__main__':
	N = 10
	eps_vec = np.linspace(0,.5,N)
	bsp_vec = np.linspace(0,10.,N)
	num_cores = multiprocessing.cpu_count()
	print(num_cores)
	results = (Parallel(n_jobs=num_cores)(delayed(processInput)(j,N,eps_vec,bsp_vec) 
		for j in range(N)))
	print(type(results))

