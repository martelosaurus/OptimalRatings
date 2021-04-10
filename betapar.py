import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# model parameters
M_up  = 1.  # upper message limit
M_dn  = 0.  # lower message limit
eps   = .1  # noise support
# sig = # standard deviation of noise (0 is Dirac; 1 is Uniform)
alpha = 10. #  
beta  = 10. #

# code parameters
N    = 1000  # #of steps in approximating step function
dx   = (M_up-M_dn)/N
tol  = 1.e-4 # error tolerance

# vector of possible m_tilde's (note: arange gives bad results)
#m_tilde = np.arange(M_dn-eps,M_up+eps,dx)
m_tilde = np.round(np.arange(M_dn-eps,M_up+eps+dx,dx),10) 

# message functions
M = np.r_[0:N+1]*dx # the step function is defined by the left point
Q_id = np.r_[0:N+1]*dx # knots for the identity message function
Q_old = np.zeros(N+1)

# guess the identity message function
Q = Q_id

# contraction
counter = 0.
while np.max(np.abs(Q-Q_old))>tol:

	# update counter
	counter += 1
	print(counter)

	# update
	Q_old = Q

	# receiver's optimal action 
	def A(mt):
		# largest and smallest possible values of q given m_tilde
		q_up = np.max(Q[M-eps<=mt])
		q_dn = np.min(Q[mt<=M+eps])
		q_diff = q_up-q_dn
		# return the conditional expection on [q_dn,q_up]
		if q_diff == 0.:
			a = q_up
		else:
			a = (q_up+q_dn)/2
			a = (stats.beta.expect(lambda x: x,(alpha,beta),
				ub=q_up,lb=q_dn,loc=q_dn,scale=q_diff,conditional=True))	
		return a #q_dn, q_up

	# vectorize A (there are many m_tilde's!)
	A = np.vectorize(A)
	Q = .5*(A(m_tilde[M_dn+eps<=m_tilde])+A(m_tilde[m_tilde<=M_up-eps]))
	Q[+0] = 0.
	Q[-1] = 1.

#print("finish update")
plt.plot(m_tilde,A(m_tilde),'-r')
plt.plot(M,Q,'-b')
plt.axis([M_dn-eps,M_up+eps,-.1,1.1])
plt.xlabel('received message $\\tilde{m}$')
plt.ylabel('action $A(\\tilde{m})$')
plt.grid()
plt.show()



