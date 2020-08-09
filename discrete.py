import numpy as np
from scipy.integrate import fixed_quad
from scipy.linalg import hankel
from functools import lru_cache

# TODO: remove this
np.set_printoptions(linewidth=160)

# -----------------------------------------------------------------------------
# Beta distribution
@np.vectorize
def _beta(e,a,b,e_bar):
	p = .5*(1.+e/e_bar)
	f = p**(a-1.)*(1.-p)**(b-1.)
	return f

@lru_cache(maxsize=None)
def _beta_wgts(a,b,e_bar):
	# NOTES: speed *should* be okay with caching, but could be a bottleneck
	return fixed_quad(lambda e : _beta(e,a,b,e_bar),-e_bar,e_bar)[0]

@np.vectorize
def beta(e,a,b,e_bar):
	"""beta distribution on [-e_bar,+e_bar]"""
	if e < -e_bar or e > e_bar:
		return 0.
	else:
		#raise Exception('Trying to evaluate error PDF out of its support')
		return _beta(e,a,b,e_bar)/_beta_wgts(a,b,e_bar)

@np.vectorize
def _alpha(i,j,M,N,a,b,e_bar,d_bar):
	y0 = 2.*d_bar*(i+M)-e_bar	# y_{i}
	y1 = y0+2.*d_bar			# y_{i+1}		
	f = lambda m_til : beta(m_til-2.*d_bar*j,a,b,e_bar)	
	return fixed_quad(f,y0,y1)[0]

# -----------------------------------------------------------------------------
# message
class Message:

	def __init__(self,M,N,a,b):
		"""
		Parameters 
		----------
		M : int
			Number of knots per message
		N : int
			Number of messages
		a, b : float
			Error parameters

		Examples
		--------

		"""
		# primitives
		self.M = M
		self.N = N

		# error distribution
		self.a = a
		self.b = b

		# knot spacing
		self.K = self.M*self.N
		self.e_bar = .5/self.N
		self.d_bar = .5/self.K

		# A matrix
		# TODO: remove I, J, A, B, C from attributes
		self.I = hankel(np.arange(-self.M,0),np.arange(-1,self.K))
		self.J = np.tile(np.arange(0,self.K+1),(self.M,1))
		self.A = _alpha(self.I,self.J,M,N,a,b,self.e_bar,self.d_bar)

		# B matrix	
		z = np.zeros(self.K)
		A_L = np.vstack((self.A[:,:-1],z))
		A_R = np.vstack((z,self.A[:,1:]))
		self.B = A_L-A_R

		# C matrix
		self.C = self.B

		# solve for breakpoints

	def alpha(self,i,j):
		return _alpha(i,j,self.M,self.N,self.a,self.b,self.e_bar,self.d_bar)
	
	def plot(self)



m2 = Message(4,2,2.,2.)
m1 = Message(4,2,1.,1.)
