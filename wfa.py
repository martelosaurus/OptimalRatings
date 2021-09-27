import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import root
from numpy.linalg import norm
from constant_f import Message

# Data (scores 0-10, no scores received 10)
m_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.009,.001])
N = len(m_wfa)
q = np.r_[0:11]

# plot distribution
plt.bar(q,100*m_wfa)
plt.xlabel("Referee score (scale of 0 to 10)")
yt = np.r_[0:30:5]
plt.yticks(yt,[str(t) + "%" for t in yt])
plt.title("Distribution of submission reviews for 2014\n"
    "Western Finance Association Meetings")
plt.savefig("wfa.pdf")

# committee reviews
m_wfa = np.hstack([0,m_wfa])
M_wfa = m_wfa.cumsum()

# calibrartion
a0 = np.hstack([1,np.zeros(N-1)])
def f(a):
    """a is a list of coefficients"""
    m = Message(N,lambda x: a.dot(np.cos(np.r_[0:N]*np.arccos(2*x-1))))
    return m.M-M_wfa

# find zero
sol = root(f,a0)
