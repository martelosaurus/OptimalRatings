from discrete import Message
from matplotlib import pyplot as plt

B = [2.,1.,0.,-.1,-.2]
while B:
	n = len(B)
	_b = B.pop()
	m = Message(M=100,N=4,b=_b)
	m.plot_msg("documentation/msg" + str(n) + ".pdf")
	m.plot_err("documentation/err" + str(n) + ".pdf",2.,6.)
