"""
Package for optimal ratings.

Organization
------------
We either assume that the optimal message function is smooth (in which case we 
compute it using the calculus of variations) or we assume that it is discrete.

The main class is 'Message'. Its children are 'Discrete' and 'Continuous'. Each
subclass has 'error' and 'plot' methods.

Model
-----
The sender and receiver agree on a message function m:[0,1]->[0,1]. The sender
privately observes q, which is uniformly distributed on [0,1]. She sends the 
receiver the message m(q). The receiver receives the message m_tilde = m(q)+e, 
where e is distributed on [-e_bar,e_bar] according to the Beta distribution. 
She then takes an action A(m_tilde). The sender and receiver incur the cost 
(q-A(m_tilde))I(q) where I:[0,1]->[0,1] is the importance function. 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import fixed_quad, dblquad
import csv

# -----------------------------------------------------------------------------
# functions that operate on lists of messages
@np.vectorize
def cost_comp(a,b):

    dc = discrete_cost(N)
    ic = identity_cost(beta,e_bar)

    return dc-ic

class Message:

	def __init__(self,func,N,I,nplot=1000):
		"""
		Parameters 
		----------
		N : int
			Number of messages
		M : 
			?
		I : callable
			Importance function
		nplot : int
			Number of knots at which to plot the message
		func : callable
			PDF of error. identity_cost assumes func has support [-e_bar,e_bar]
		"""
		self.N = N
		self.I = I
		self.M = D(I,self.N) 
		self.nplot = nplot

def comp_plot(m1,m2,a_max=2,b_max=2,n_plot=40):
	"""
	Message comparison plot

	Parameters
	----------
	m1, m2 : Message
		Messages to compare
	a_max, b_max : float
		Maximum 'a' and 'b' parameters (the minimums are zero for both)
	n_plot : int
		Square root of the number of (a,b) knots at which to compare the 
		message errors
	"""

	# plotting vectors, matrices
	a_vec = np.linspace(0.,a_max,n_plot)
	b_vec = np.linspace(0.,b_max,n_plot)
	A, B = np.meshgrid(a_vec,b_vec)
	#I = cost_comp(A,B)
	@np.vectorize
	def heaviside(x):
		if x>0:
			return x
		else:
			return np.nan
	I = pd.read_csv('Imat.csv',header=None).to_numpy()
	#I = heaviside(I)

	# main plot calls
	fig, axs = plt.subplots()
	axs.contourf(A,B,I,[0.,1.],colors=['lightgrey'])
	axs.plot(a_vec,a_vec,color='grey',linestyle='-',linewidth=.5)
	axs.plot(1.,1.,'ok')
	axs.annotate('Identity $\\succsim$ Discrete',(1.2,1.5))
	axs.annotate('Discrete $\\succsim$ Identity',(.2,.5))
	axs.annotate('Uniform Distribution',(1.05,.9))

	# axes
	axs.set_xticks([0.,1.,a_max])
	axs.set_yticks([0.,1.,b_max])
	axs.set_xticklabels(['$0$','$1$','$2$'])
	axs.set_yticklabels(['$0$','$1$','$2$'])
	axs.set_xlabel('$\\alpha$')
	axs.set_ylabel('$\\beta$')
	axs.set_aspect('equal','box')
	axs.grid(color='grey',linestyle='-',linewidth=.5)

	# figure
	fig.tight_layout()
	fig.show()


def mplot(M,nplot=1000):
	"""
	Plot a list of message functions
	
	Parameters 
	----------
	M : list (of Messages)
		List of messages
	nplot : int
		Number of knots at which to plot the message
	"""
	I = M[0].I
	q = np.linspace(0.,1.,nplot)
	cplt = np.zeros(nplot)
	for n in range(0,nplot):
		cplt[n] = _cts_msg_fun(I,q[n])
	plt.plot(q,cplt,'-r')
	for m in M:
		plt.plot(q,(np.digitize(q,m.M)-1.)/(m.N-1.),'-b')
	plt.plot(q,q,'--k')
	plt.plot(q,cplt,'-r')
	plt.axis([0.,1.,0.,1.])
	plt.legend(['asymptotic message','optimal discrete'])
	plt.xlabel('state $(q)$')
	plt.ylabel('message $(m)$')
	plt.grid()
	plt.show()

def mdrop(M,fname='messages.csv',ndrop=1000):
	"""
	Writes a list of message function to disk
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]

	"""
	I = M[0].I
	q = np.linspace(0.,1.,ndrop)
	cplt = np.zeros(ndrop)
	with open(fname,'w') as csvfile:
		writer = csv.writer(csvfile)
		R = ['q','I','m_inf']
		for m in M:
			R = R + ['m_'+str(m.N)]
		writer.writerow(R)
		for n in range(0,ndrop):
			R = [q[n],I(q[n]),_cts_msg_fun(I,q[n])]
			for m in M:
				R = R + [np.digitize(q[n],m.M)/(1.*m.N)] 
			writer.writerow(R)

def eplot(m,nplot=4):
	"""
	Plots the error against the N or e_bar
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]
	"""
	E = np.zeros(nplot)
	for n in range(0,nplot):
		E[n] = Message(10**(n+1),I).error()
	plt.plot(np.r_[1:(nplot+1)],np.log10(E))
	plt.xlabel('$\\log_{10}$(Number of Messages)')
	plt.ylabel('$\\log_{10}$(Relative $L^{2}$ Error)')
	plt.grid()
	plt.show()
