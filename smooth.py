class Smooth(Message):

	def __init__(self):

	def _cts_msg_fun(I,q):
		"""
		optimal [C]ontinuous message function
		
		Parameters
		----------
		I : callable
			Importance function I:[0,1]->[0,1]

		"""
		return _Q(I,0.,(1./3),0.,q)/_Q(I,0.,(1./3),0.,1.)

	def cost(self):
		"""Returns the L^{2} error"""
		E = 0.
		I0 = integrate.quad(lambda x: _cts_msg_fun(self.I,x)**2,0.,1.)[0]
		for n in range(0,self.N):
			f = lambda x: (n/(1.*self.N)-_cts_msg_fun(self.I,x))**2
			E = E + integrate.quad(f,self.M[n],self.M[n+1])[0]
		return np.sqrt(E/I0)

		
