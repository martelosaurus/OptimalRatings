from discrete import Message
from matplotlib import pyplot as plt

B = [-.2,-.1,0.,1.,2.,4.]
while B:
	n = len(B)
	_b = B.pop()
	m = Message(M=100,N=4,b=_b)
	m.plot_msg("../latex/Model_Figures/msg" + str(n) + ".pdf",title=False)
	m.plot_err("../latex/Model_Figures/err" + str(n) + ".pdf",2.,6.,title=False)
