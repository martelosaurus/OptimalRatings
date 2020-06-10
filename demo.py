"""This is a sandbox script for the discrete message module"""
import numpy as np
import pandas as pd
from optrat import *
from scipy.integrate import fixed_quad
from functools import lru_cache
#import dm as dm


# beta test
N = 10
e_bar = .5/N

def uniform_identity_cost():
    return (1./3.)*(e_bar**2.)*(1.-e_bar)

def uniform_discrete_cost():
    return (1./3.)*(e_bar**2.)/(1.+2.*e_bar)**2.

@np.vectorize
def cost_comp(a,b):

    @np.vectorize
    def _beta(e):
        p = .5*(1.+e/e_bar)
        f = p**(a-1.)*(1.-p)**(b-1.)
        return f

    w = fixed_quad(_beta,-e_bar,e_bar)[0]

    @np.vectorize
    def beta(e):
        if e < -e_bar or e > e_bar:
            raise Exception('Trying to evaluate error PDF out of its support')
        return _beta(e)/w

    dc = discrete_cost(N)
    ic = identity_cost(beta,e_bar)

    #print('{:.5f}, {:.5f}'.format(dc,ic))
    return dc-ic

# -----------------------------------------------------------------------------
# plotting
n_plot = 40 
a_max = 2.
b_max = 2.

# plotting vectors, matrices
a_vec = np.linspace(0.,a_max,n_plot)
b_vec = np.linspace(0.,b_max,n_plot)
A, B = np.meshgrid(a_vec,b_vec)
#I = cost_comp(A,B)
@np.vectorize
def heaviside(x):
    if x>0:
        return x
    else:
        return np.nan
I = pd.read_csv('Imat.csv',header=None).to_numpy()
#I = heaviside(I)

# main plot calls
fig, axs = plt.subplots()
axs.contourf(A,B,I,[0.,1.],colors=['lightgrey'])
axs.plot(a_vec,a_vec,color='grey',linestyle='-',linewidth=.5)
axs.plot(1.,1.,'ok')
axs.annotate('Identity $\\succsim$ Discrete',(1.2,1.5))
axs.annotate('Discrete $\\succsim$ Identity',(.2,.5))
axs.annotate('Uniform Distribution',(1.05,.9))

# axes
axs.set_xticks([0.,1.,a_max])
axs.set_yticks([0.,1.,b_max])
axs.set_xticklabels(['$0$','$1$','$2$'])
axs.set_yticklabels(['$0$','$1$','$2$'])
axs.set_xlabel('$\\alpha$')
axs.set_ylabel('$\\beta$')
axs.set_aspect('equal','box')
axs.grid(color='grey',linestyle='-',linewidth=.5)

# figure
fig.tight_layout()
fig.show()
