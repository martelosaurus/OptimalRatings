import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import fixed_quad
from scipy.optimize import root
from scipy.special import erf
from numpy.linalg import norm

# parameters
s = 1. # shape parameters

# Data (scores 0-10, no scores received 10)
x_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.009,.001])
N = len(x_wfa)
q = np.r_[0:11]

# plot distribution
plt.bar(q,100*x_wfa)
plt.xlabel("Referee score (scale of 0 to 10)")
yt = np.r_[0:30:5]
plt.yticks(yt,[str(t) + "%" for t in yt])
plt.title("Distribution of submission reviews for 2014\n"
    "Western Finance Association Meetings")
plt.savefig("wfa.pdf")
plt.close()

# committee reviews
x_wfa = np.hstack([0,x_wfa])
X_wfa = x_wfa.cumsum()

# calibrartion
a0 = np.ones(N+1)/N
xk = np.linspace(0,1,N+1)

def _Q(i,j,a,b,A=.1,_n=11):
    return fixed_quad(lambda q: (q**i)*np.exp(-A*(q-xk[j])**2),a,b,n=_n)[0]
