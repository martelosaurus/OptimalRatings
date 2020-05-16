"""
This is a sandbox script for the discrete message module
"""
import numpy as np
from optrat import *
from scipy.integrate import fixed_quad
from functools import lru_cache
#import dm as dm

if False:
    # importance function(s)
    I = {
            'i1' : lambda x: x**3., 
            'i2' : lambda x: x**(-1.5),
            'i3' : lambda x: (6.*(x-.5)**2.+.5)**3.
        }

    # for each importance function, plot discrete messages of size 5, 20, 100 
    for i in I:
        M = []
        for n in [5,20,100]:
            M = M + [dm.Message(n,I[i])]
        dm.mplot(M)

# beta test
N = 10
e_bar = .5/N

@np.vectorize
def cost_comp(a,b):

    print('({a:.3f}, {b:.3f})'.format(a=a,b=b))
    @np.vectorize
    def _beta(e):
        p = .5*(1.+e/e_bar)
        f = p**(a-1.)*(1.-p)**(b-1.)
        return f

    w = fixed_quad(_beta,-e_bar,e_bar)[0]

    @np.vectorize
    def beta(e):
        if e < -e_bar or e > e_bar:
            raise Exception('Trying to evaluate error PDF outside of its support')
        return _beta(e)/w

    return discrete_cost(beta,N)-identity_cost(beta,e_bar)

# -----------------------------------------------------------------------------
# plotting

n_plot = 100 
a_max = 2.
b_max = 2.

# plotting vectors, matrices
a_vec = np.linspace(1.,a_max,n_plot)
b_vec = np.linspace(1.,b_max,n_plot)
A, B = np.meshgrid(a_vec,b_vec)
I = cost_comp(A,B)

# main plot calls
fig, axs = plt.subplots()
axs.contour(A,B,I)
axs.plot(a_vec,a_vec,'k--')
axs.plot(1.,1.,'ok')

# axes
axs.set_xticks([0.,1.])
axs.set_yticks([0.,1.])
axs.set_xticklabels(['$0$','$1$'])
axs.set_yticklabels(['$0$','$1$'])
axs.set_xlabel('$\\alpha$')
axs.set_ylabel('$\\beta$')
axs.set_aspect('equal','box')

# figure
fig.tight_layout()
fig.show()
