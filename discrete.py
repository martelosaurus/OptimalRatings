import numpy as np
from scipy.integrate import fixed_quad
from scipy.linalg import hankel
from functools import lru_cache

# -----------------------------------------------------------------------------
# Beta distribution
@np.vectorize
def _beta(e,a,b,e_bar):
    p = .5*(1.+e/e_bar)
    f = p**(a-1.)*(1.-p)**(b-1.)
    return f

@lru_cache(maxsize=None)
def _beta_wgts(a,b,e_bar):
    # NOTES: speed *should* be okay with caching, but could be a bottleneck
    return fixed_quad(lambda e : _beta(e,a,b,e_bar),-e_bar,e_bar)[0]

@np.vectorize
def beta(e,a,b,e_bar):
    """beta distribution on [-e_bar,+e_bar]"""
    if e < -e_bar or e > e_bar:
        raise Exception('Trying to evaluate error PDF out of its support')
    return _beta(e,a,b,e_bar)/_beta_wgts(a,b,e_bar)

@np.vectorize
def _alpha(i,j,a,b,e_bar):
    return fixed_quad(lambda e : beta(e,a,b,e_bar),-e_bar,e_bar)[0]

# -----------------------------------------------------------------------------
# message
class Message:

    def __init__(self,M,N,a,b):
        """
        Parameters 
        ----------
        M : int
            Number of knots per message
        N : int
            Number of messages
        a, b : float
            Error parameters

        Examples
        --------

        # primitives
        """
        self.M = M
        self.N = N

        # knot spacing
        self.K = self.M*self.N
        self.e_bar = .5/self.N
        self.d_bar = .5/self.K

        # index matrices
        # TODO: remove I, J, A, B, C from attributes
        self.I = np.tile(np.arange(0,self.K+1),(self.M,1))
        self.J = hankel(np.arange(-self.M,0),np.arange(-1,self.K))
        self.A = _alpha(self.I,self.J,self.a,self.b,self.e_bar)
        self.B = None
        self.C = None
        
        # solve

