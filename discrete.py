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

		# index matrices
		# TODO: remove I, J, A, B, C from attributes
		self.I = np.tile(np.arange(0,self.K+1),(self.M,1))
		self.J = hankel(np.arange(-self.M,0),np.arange(-1,self.K))
		self.A = _alpha(self.I,self.J,M,N,a,b,self.e_bar,self.d_bar)
		self.B = None
		self.C = None

		# check
		ii = np.arange(-M,self.K)
		jj = np.arange(0,self.K)
		II, JJ = np.meshgrid(ii,jj)
		self.II = II
		self.JJ = JJ
		self.AA = _alpha(II,JJ,M,N,a,b,self.e_bar,self.d_bar)
		
		X = list(zip(self.I.flatten(),self.J.flatten()))
		@np.vectorize
		def _alpha_check(i,j):
			if (i,j) in X:
				return 1.
			else:
				return 0.
		self.XX = _alpha_check(II,JJ)
		
		# solve

