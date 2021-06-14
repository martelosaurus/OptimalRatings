
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
# matplotlib imports
from matplotlib import rc
from matplotlib import pyplot as plt

# numpy imports
import numpy as np

# scipy imports
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
    plt.plot(q,cplt,color="tab:orange")
    for m in M:
        plt.plot(q,(np.digitize(q,m.M)-1.)/(m.N-1.),color="tab:blue")
    plt.plot(q,q,'--k')
    plt.plot(q,cplt,color="tab:orange")
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
