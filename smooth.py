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

