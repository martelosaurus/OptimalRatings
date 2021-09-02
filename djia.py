import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.neighbors import KernelDensity

# load data
m = pd.read_csv("../data/HistoricalPrices.csv")

# instantiate and fit KDE model
kde = KernelDensity(bandwidth=1.0, kernel='gaussian')
p = m[" Close"]
r = p.to_numpy()/p.shift().to_numpy()-1.
r = r[1:]
r.sort()
f = kde.fit(r[:,None])
r_plt = np.linspace(-.1,.1,100)
logprob = kde.score_samples(r_plt[:,None])
plt.fill_between(r_plt, np.exp(logprob), alpha=0.5)
plt.plot(r, np.full_like(r, -0.01), '|k', markeredgewidth=1)
plt.show()

if False:

	# score_samples returns the log of the probability density
	logprob = kde.score_samples(x_d[:, None])

	plt.fill_between(x_d, np.exp(logprob), alpha=0.5)
	plt.plot(x, np.full_like(x, -0.01), '|k', markeredgewidth=1)
	plt.ylim(-0.02, 0.22)

# candidate message function

