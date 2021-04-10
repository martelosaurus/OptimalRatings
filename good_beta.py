import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, integrate

# driver for identity messages on the unit interval
def int_identity(eps,a):
	# middle
	I_m = (1-2*eps)*stats.beta.expect(lambda e : e**2,(a,a),-eps,2*eps) 
	# bottom-left
	i_b = lambda e, q: ((q-(q+e+eps+0.)/2)**2)*stats.beta.pdf(e,a,a,-eps,2*eps)
	I_b = integrate.dblquad(i_b,0,2*eps,lambda q: -eps,lambda q: eps-q)[0]
	# top-right
	i_t = lambda e, q: ((q-(q+e-eps+1.)/2)**2)*stats.beta.pdf(e,a,a,-eps,2*eps)
	I_t = integrate.dblquad(i_t,1-2*eps,1.,lambda q: 1.-eps-q,lambda q: eps)[0]
	return I_m+I_b+I_t
int_identity = np.vectorize(int_identity)

# driver for discrete messages on the unit interval
def int_discrete(eps,a):
	return ((eps/(2*eps+1))**2)/3
int_discrete = np.vectorize(int_discrete)

# parameter vector
eps = .2
n   = 100
f = lambda e, a: stats.beta.pdf(e,a,a,-eps,2*eps)
N = 30
bsp_vec = np.linspace(1.,10.,N)
cost_discrete = int_discrete(eps*np.ones(N),bsp_vec)
cost_identity = int_identity(eps*np.ones(N),bsp_vec)

# PLOT 1: beta distributions
e   = np.linspace(-eps,eps,n)
with open('beta_pdfs.csv','w') as csvfile:
	contourwriter = csv.writer(csvfile)
	contourwriter.writerow(['e','PDF_a1','PDF_a2','PDF_a10'])
	for j in range(0,n-1):
		row = np.hstack([e[j],f(e[j],1.),f(e[j],2.),f(e[j],10.)])
		contourwriter.writerow(row)

# PLOT 1: discrete vs identity
with open('discrete_identity.csv','w') as csvfile:
	contourwriter = csv.writer(csvfile)
	contourwriter.writerow(['a','cost_discrete','cost_identity'])
	for j in range(0,N-1):
		row = np.hstack([bsp_vec[j],cost_discrete[j],cost_identity[j]])
		contourwriter.writerow(row)