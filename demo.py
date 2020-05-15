"""
This is a sandbox script for the discrete message module
"""
import numpy as np
from optrat import *
from scipy.integrate import fixed_quad
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
a = 2.
b = 2.
@np.vectorize
def beta(x):
    p = (2.*x-1.)*e_bar 
    f = p**(a-1.)*(1.-p)**(b-1.)
    w = fixed_quad(f,-e_bar,e_bar)
    return f/w

# -----------------------------------------------------------------------------
# plotting
if False:

    # plotting vectors, matrices
    a_vec = np.linspace(0.,a_max,n_plot)
    b_vec = np.linspace(0.,b_max,n_plot)
    A, B = np.meshgrid(a_vec,b_vec)

    # main plot calls
    fig, axs = plt.subplots()
    axs.contour(A,B,A**2.+B**2.,[1.])
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
        
