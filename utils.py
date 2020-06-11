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

def _Q(I,i,j,a,b,_n=11):
	"""
	[Q]uadrature: returns \int_{a}^{b}{q^{i}I(q)^{j}dq)
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]
	_n : int
		Number of quadrature points

	"""
	return fixed_quad(lambda q: (q**i)*(I(q)**j),a,b,n=_n)[0]

def _A(I,x1,x2):
	"""
	Optimal [A]ction given that x1 < q < x2
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]
	x1, x2 : float, float
		Given x1 < q < x2
	"""
	if x1 == x2:
		return x1
	else:
		return _Q(I,1.,1.,x1,x2)/_Q(I,0.,1.,x1,x2)
	
def _X(I,x1,x2):
	"""
	ne[X]t step
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]

	"""
	if x1 == x2:
		return x1
	else:
		f = lambda x3: x2-(_A(I,x1,x2)+_A(I,x2,x3))/2.
		return optimize.newton(f,x2)

def _S(I,x1,N):
	"""
 	non-linear [S]hooting
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]

	"""
	x = np.zeros(N+1)
	x[1] = x1
	for n in range(0,N-1): 
		x[n+2] = _X(I,x[n],x[n+1]) 
	return x

def D(I,N):
	"""
	optimal [D]iscrete message function
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]

	"""
	x1s = optimize.newton(lambda x1: _S(I,x1,N)[-1]-1.,1./N)
	return _S(I,x1s,N) 
