class Discrete(Message):

	def __init__(self):

		#self.knots = _knots(x)

	def cost(self,N):
		"""
		Cost of discrete message function 
		
		Parameters
		----------
		N : int
			There are N+1 messages

		Notes
		-----
		It's important to note that the discrete message function has the same
		cost regardless of the error distribution. Morever, if the importance 
		function is constant, then the cost is given by

			cost[discrete] = (1./3.)*(e_bar**2.)/(1.+2.*e_bar)**2.

		where e_bar = 1./(2.*N). 
		"""
		return 1. 
