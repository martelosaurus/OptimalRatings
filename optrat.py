"""
Package for optimal ratings.

Organization
------------
We either assume that the optimal message function is smooth (in which case we 
compute it using the calculus of variations) or we assume that it is discrete.

The main class is 'Message'. Its children are 'Discrete' and 'Continuous'. Each
subclass has 'error' and 'plot' methods.

Model
-----
The sender and receiver agree on a message function m:[0,1]->[0,1]. The sender
privately observes q, which is uniformly distributed on [0,1]. She sends the 
receiver the message m(q). The receiver receives the message m_tilde = m(q)+e, 
where e is distributed on [-e_bar,e_bar] according to the Beta distribution. 
She then takes an action A(m_tilde). The sender and receiver incur the cost 
(q-A(m_tilde))I(q) where I:[0,1]->[0,1] is the importance function. 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import fixed_quad, dblquad
import csv

# -----------------------------------------------------------------------------
# functions that operate on lists of messages
@np.vectorize
def cost_comp(a,b):

    dc = discrete_cost(N)
    ic = identity_cost(beta,e_bar)

    return dc-ic

class Message:

    def __init__(self,func,N,I):
        """
        Parameters 
        ----------
        N : int
            Number of messages
        M : int
            Number of knots per message
        I : callable
            Importance function
        func : callable
            PDF of error. identity_cost assumes func has support [-e_bar,e_bar]

        Examples
        --------

        """
        self.N = N
        self.I = I
        self.M = D(I,self.N) 

        self.K = self.M*self.N
        self.e_bar = .5/self.N
        self.d_bar = .5/self.K

