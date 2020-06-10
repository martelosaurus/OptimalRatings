class Discrete(Message):

	def __init__(self):

	def discrete_cost(N):
		"""
		Cost of discrete message function with constant importance function
		
		Parameters
		----------
		N : int
			There are N+1 messages

		Notes
		-----
		The discrete message function has the same cost regardless of the error
		distribution. 
		"""
		e_bar = 1./(2.*N) 
		return (1./3.)*(e_bar**2.)/(1.+2.*e_bar)**2.
