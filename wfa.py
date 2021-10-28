import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import fixed_quad
from scipy.optimize import root
from scipy.sparse import spdiags
from numpy.linalg import norm

# Data (scores 0-10, no scores received 10)
x_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.01,.0])
#x_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.005,.005])
N = len(x_wfa)
q = np.r_[0:11]

# plot distribution
plt.bar(q,100*x_wfa)
plt.xlabel("Referee score (scale of 0 to 10)")
xt = np.r_[0:11]
yt = np.r_[0:30:5]
plt.xticks(xt,xt)
plt.yticks(yt,[str(t) + "%" for t in yt])
plt.title("Distribution of submission reviews for 2014\n"
    "Western Finance Association Meetings")
#plt.show()
plt.savefig("wfa_hist.pdf")
plt.close()

# committee reviews
X_wfa = x_wfa.cumsum()

# halfpoints
g_wfa_x = .5*(X_wfa[1:]+X_wfa[:-1])
g_wfa_y = (1/(X_wfa[1:]-X_wfa[:-1]))**3

# plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot((1+q)/N,q,"--",color="tab:blue",linewidth=2)
ax1.plot(X_wfa,q,"-",color="tab:blue",linewidth=2)
ax2.plot(g_wfa_x,g_wfa_y,'.k')
ax1.set_xlabel("Paper quality percentile")
ax1.set_ylabel("Referee score")
ax2.set_ylabel("Density")
xt = np.linspace(0,1,N)
yt = np.r_[0:11]
ax1.set_xticks(xt)
ax1.set_xticklabels([str(t)+"%" for t in np.r_[0:110:10]])
ax1.set_yticks(yt)
ax2.set_yticklabels(yt)
ax1.legend(["Straight talk","WFA reviews"])
ax2.legend(["Inferred importance"],loc="upper center")
plt.title("WFA Program Committee Reviews")
plt.savefig("wfa_dist.pdf")
