import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate

n = 1000
act = np.zeros(n)
sig = np.linspace(.1,10.,n)
for j in range(1,n):
	act[j] = integrate.quad(lambda x: np.exp(-(x/sig[j])**2)*x**2,-1,1.)[0]
plt.plot(sig,act)
plt.grid()
plt.show()