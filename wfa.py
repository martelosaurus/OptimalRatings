import numpy as np
from matplotlib import pyplot as plt

# ------------------------------------------------------------------------------
# Data 
# (scores 0-10, no scores received 10)
x_wfa = np.array([.05,.20,.23,.17,.14,.10,.06,.02,.02,.01,.0])
X_wfa = x_wfa.cumsum()
N = len(x_wfa)
m = np.r_[0:N]

# ------------------------------------------------------------------------------
# message distribution
plt.bar(m,100*x_wfa)
# x-axis
plt.xlabel("Referee score (scale of 0 to 10)")
m = np.r_[0:N]
plt.xticks(m,m)
# y-axis
yt = np.r_[0:30:5]
plt.yticks(yt,[str(t) + "%" for t in yt])
plt.title("Distribution of submission reviews for 2014\n"
    "Western Finance Association Meetings")
# save
plt.savefig("../text/Model_Figures/wfa_hist.pdf")
plt.close()

# ------------------------------------------------------------------------------
# quality distribution (at halfpoints)
g_wfa_x = .5*(X_wfa[1:]+X_wfa[:-1])
g_wfa_y = (.1/(X_wfa[1:]-X_wfa[:-1]))**3

# plot
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot((1+m)/N,m,"--",color="tab:blue",linewidth=2)
ax1.plot(X_wfa,m,".-",color="tab:blue",linewidth=2)
ax2.plot(g_wfa_x,g_wfa_y,'.-k')
ax1.set_xlabel("Paper quality percentile ($q$)")
ax1.set_ylabel("Referee score")
ax2.set_ylabel("Importance")
xt = np.linspace(0,1,N)
yt = np.r_[0:N]
ax1.set_xticks(xt)
ax1.set_xticklabels([str(t)+"%" for t in np.r_[0:110:10]])
ax1.set_yticks(yt)
ax2.set_yticklabels(yt)
ax1.axis([0,1,-1,11])
ax2.axis([0,1,-1,11])
ax1.legend(["Straight talk","WFA reviews"])
ax2.legend(["Inferred importance"],loc="upper center")
plt.title("WFA Program Committee Reviews")
plt.savefig("../text/Model_Figures/wfa_dist.pdf")
