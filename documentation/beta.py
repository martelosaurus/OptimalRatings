import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import uniform
from error import error
from scipy.integrate import fixed_quad

# -----------------------------------------------------------------------------
# parameters
n_plot = 100 
a_max = 2.
b_max = 2.
e_bar = .1 # raise an exception if this doesnt divide 1

# -----------------------------------------------------------------------------
# functions
@np.vectorize
def q_dn(m_til):
    return max(m_til-e_bar,0.)

@np.vectorize
def q_up(m_til):
    return min(m_til+e_bar,1.)

@np.vectorize
def identity_cost(a,b):
    """Cost of identity message function"""

    Q = uniform()
    E = error(a,b,e_bar) 

    @np.vectorize
    def A(m_til): # optimal action given a received message m_til
        return E.expect(lambda e: e,m_til-q_dn(m_til),m_til+q_up(m_til))

    @np.vectorize
    def I(q):
        return E.expect(lambda e: (q-A(q+e))**2.,-e_bar,e_bar)

    # TODO: joint expectation
    return Q.expect(I,0.,1.)

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
