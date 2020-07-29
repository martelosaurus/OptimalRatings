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
def beta(e,a,b):
    """beta distribution on [-e_bar,+e_bar]"""
    if e < -e_bar or e > e_bar:
        raise Exception('Trying to evaluate error PDF out of its support')
    return _beta(e,a,b,e_bar)/_beta_wgts(a,b,e_bar)

# -----------------------------------------------------------------------------
# quadrature
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
