import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import fixed_quad
from scipy.optimize import root
from scipy.sparse import spdiags
from numpy.linalg import norm

# parameters
s = 1. # shape parameters

# Data (scores 0-10, no scores received 10)
x_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.01,.0])
N = len(x_wfa)
q = np.r_[0:11]

# plot distribution
plt.bar(q,100*x_wfa)
plt.xlabel("Referee score (scale of 0 to 10)")
yt = np.r_[0:30:5]
plt.yticks(yt,[str(t) + "%" for t in yt])
plt.title("Distribution of submission reviews for 2014\n"
    "Western Finance Association Meetings")
plt.show()
plt.savefig("wfa.pdf")
plt.close()

# committee reviews
x_wfa = np.hstack([0,x_wfa])
X_wfa = x_wfa.cumsum()
g_wfa_x = .5*(X_wfa[1:]+X_wfa[:-1])
g_wfa_y = 
plt.xlabel("$q$")
plt.ylabel("$m(q)$")
plt.plot(X_wfa,X_wfa)
plt.plot(X_wfa,x_plt)
plt.grid()
plt.show()
