# matplotlib imports
from matplotlib import rc
from matplotlib import pyplot as plt

# numpy imports
import numpy as np

# scipy imports
from scipy import integrate, optimize
from scipy.integrate import fixed_quad, dblquad
import csv

def _fixed_quad(*args,**kwargs):
    """Wrapper for fixed_quad that returns the integral (without the error)"""
    return fixed_quad(*args,**kwargs)[0]

def fixed_dblquad_tri(func,a,b,gfun,hfun):
    """
    Homebrewed double quadrature using scipy.integrate.fixed_quad.
    Parameters 
    ----------
    See help(scipy.integrate.dblquad)
    """
    f = lambda x: _fixed_quad(lambda y: func(y,x),gfun(x),hfun(x))
    return _fixed_quad(np.vectorize(f),a,b)

def identity_cost(func,e_bar,tol=1.e-10):
    """
    Cost of identity message function with constant importance function
	'identity_cost' stands alone because it needs to run faster than 
	Smooth.cost. 
    
    Parameters
    ----------
    func : callable
        PDF of error. identity_cost assumes func has support [-e_bar,e_bar]
    e_bar : float
        Support parameter for the error distribution (see above).
    tol : float
        Error PDF should integrate to unity (+/- tol)
    Notes
    -----
    The action function kinks at 1-e_bar and 1+e_bar. Therefore, identity_cost 
    breaks the integration problem into three parts. The action function is 
    computed using scipy.integrate.fixed_quad, which caches weights/knots, so
    should be fast. 
    There are 9 anonymous functions. Guido would not approve.
    """

    if np.abs(_fixed_quad(func,-e_bar,e_bar)-1.)>tol:
        raise Exception('PDF of error does not integrate to one')

    # quasi-expectation
    @np.vectorize
    def A(m_til,a,b):
        if a == b:
            return a
        else:
            I0 = _fixed_quad(lambda x: func(m_til-x)*x**0.,a,b)
            I1 = _fixed_quad(lambda x: func(m_til-x)*x**1.,a,b)
            return I1/I0

    # non-constant m_tilde limits of integration (see dblquad documentation)
    gfun = lambda m_til : m_til-e_bar # lower
    hfun = lambda m_til : m_til+e_bar # upper 

    # integration over [1-e_bar,1+e_bar]
    a_up = lambda m_til: A(m_til,m_til-e_bar,1.)
    f_up = lambda q, m_til: func(m_til-q)*(q-a_up(m_til))**2.
    z_up = fixed_dblquad_tri(f_up,1.-e_bar,1.+e_bar,gfun,lambda _:1.)

    # integration over [+e_bar,1-e_bar]
    a_md = lambda m_til: A(m_til,m_til-e_bar,m_til+e_bar)
    f_md = lambda q, m_til: func(m_til-q)*(q-a_md(m_til))**2.
    z_md = fixed_dblquad_tri(f_md,e_bar,1.-e_bar,gfun,hfun)

    # integration over [-e_bar,+e_bar]
    a_dn = lambda m_til: A(m_til,0.,m_til+e_bar)
    f_dn = lambda q, m_til: func(m_til-q)*(q-a_dn(m_til))**2.
    z_dn = fixed_dblquad_tri(f_dn,-e_bar,e_bar,lambda _:0.,hfun)

    return z_up + z_md + z_dn

def discrete_cost(N):
    """
    Cost of discrete message function with constant importance function
    
    Parameters
    ----------
    N : int
        There are N+1 messages
	Notes
	-----
	The discrete message function has the same cost regardless of the error 
	distribution. 
    """
    e_bar = 1./(2.*N) 
    return (1./3.)*(e_bar**2.)/(1.+2.*e_bar)**2.
