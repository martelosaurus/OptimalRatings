import time
import numpy as np
from scipy.integrate import fixed_quad, dblquad

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
    #return _fixed_quad(np.vectorize(f),a,b)

# demo
gfun = lambda x: 0.
hfun = lambda x: 1.-x
a = 0.
b = 1.
@np.vectorize
def func(y,x):
    return x**2.+y**2.
t1 = time.time()
for j in range(0,1000):
    Q1 = fixed_dblquad_tri(func,a,b,gfun,hfun)
t2 = time.time()
print(t2-t1)
t1 = time.time()
for j in range(0,1000):
    Q2 = dblquad(func,a,b,gfun,hfun)[0]
t2 = time.time()
print(t2-t1)
