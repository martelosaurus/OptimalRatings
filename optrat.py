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
where e is distributed on [-e_bar,e_bar] according to the PDF f. She then takes
an action A(m_tilde). The sender and receiver incur the cost (q-A(m_tilde))I(q)
where I:[0,1]->[0,1] is the importance function. 
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import fixed_quad, dblquad
import csv

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

#------------------------------------------------------------
class Message:

    def __init__(self):

class Continuous(Message):

    def __init__(self):

class Discrete(Message):

    def __init__(self):

class Identity(Message):

    def __init__(self):

#------------------------------------------------------------

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

def _cts_msg_fun(I,q):
    """
    optimal [C]ontinuous message function
    
    Parameters
    ----------
    I : callable
            Importance function I:[0,1]->[0,1]

    """
    return _Q(I,0.,(1./3),0.,q)/_Q(I,0.,(1./3),0.,1.)

class Message():

    def __init__(self,N,I,nplot=1000):
        """
        Parameters 
        ----------
        N : int
                ?
        M : 
                ?
        I : callable
                Importance function
        nplot : int
                Number of knots at which to plot the message
        """
        self.N = N
        self.I = I
        self.M = D(I,self.N) 
        self.nplot = nplot

    def error(self):
        """Returns the L^{2} error"""
        E = 0.
        I0 = integrate.quad(lambda x: _cts_msg_fun(self.I,x)**2,0.,1.)[0]
        for n in range(0,self.N):
            f = lambda x: (n/(1.*self.N)-_cts_msg_fun(self.I,x))**2
            E = E + integrate.quad(f,self.M[n],self.M[n+1])[0]
        return np.sqrt(E/I0)

def mplot(M,nplot=1000):
    """
    Plot a list of message functions
    
    Parameters 
    ----------
    M : list (of Messages)
            List of messages
    nplot : int
            Number of knots at which to plot the message
    """
    I = M[0].I
    q = np.linspace(0.,1.,nplot)
    cplt = np.zeros(nplot)
    for n in range(0,nplot):
            cplt[n] = _cts_msg_fun(I,q[n])
    plt.plot(q,cplt,'-r')
    for m in M:
            plt.plot(q,(np.digitize(q,m.M)-1.)/(m.N-1.),'-b')
    plt.plot(q,q,'--k')
    plt.plot(q,cplt,'-r')
    plt.axis([0.,1.,0.,1.])
    plt.legend(['asymptotic message','optimal discrete'])
    plt.xlabel('state $(q)$')
    plt.ylabel('message $(m)$')
    plt.grid()
    plt.show()

def mdrop(M,fname='messages.csv',ndrop=1000):
	"""
	Writes a list of message function to disk
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]

	"""
	I = M[0].I
	q = np.linspace(0.,1.,ndrop)
	cplt = np.zeros(ndrop)
	with open(fname,'w') as csvfile:
		writer = csv.writer(csvfile)
		R = ['q','I','m_inf']
		for m in M:
			R = R + ['m_'+str(m.N)]
		writer.writerow(R)
		for n in range(0,ndrop):
			R = [q[n],I(q[n]),_cts_msg_fun(I,q[n])]
			for m in M:
				R = R + [np.digitize(q[n],m.M)/(1.*m.N)] 
			writer.writerow(R)

def eplot(I,nplot=4):
	"""
	Plots the error against the number of messages
	
	Parameters
	----------
	I : callable
		Importance function I:[0,1]->[0,1]
	"""
	E = np.zeros(nplot)
	for n in range(0,nplot):
		E[n] = Message(10**(n+1),I).error()
	plt.plot(np.r_[1:(nplot+1)],np.log10(E))
	plt.xlabel('$\\log_{10}$(Number of Messages)')
	plt.ylabel('$\\log_{10}$(Relative $L^{2}$ Error)')
	plt.grid()
	plt.show()

# -----------------------------------------------------------------------------
# functions
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