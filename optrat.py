import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize
from scipy.integrate import fixed_quad, dblquad
import csv

if False:
    def Q(I,m,n,a,b):
        """[Q]uadrature"""
        return integrate.fixed_quad(lambda q: (q**m)*(I(q)**n),a,b,n=11)[0]

    def A(I,x1,x2):
        """[A]ction function"""
        if x1 == x2:
            return x1
        else:
            return Q(I,1.,1.,x1,x2)/Q(I,0.,1.,x1,x2)
        
    def X(I,x1,x2):
        """ne[X]t step"""
        if x1 == x2:
            return x1
        else:
            f = lambda x3: x2-(A(I,x1,x2)+A(I,x2,x3))/2.
            return optimize.newton(f,x2)

    def S(I,x1,N):
        """ non-linear [S]hooting"""
        x = np.zeros(N+1)
        x[1] = x1
        for n in range(0,N-1): 
            x[n+2] = X(I,x[n],x[n+1]) 
        return x

    def D(I,N):
        """optimal [D]iscrete message function"""
        x1s = optimize.newton(lambda x1: S(I,x1,N)[-1]-1.,1./N)
        print(S(I,x1s,N))
        return S(I,x1s,N) 

    def C(I,q):
        """optimal [C]ontinuous message function"""
        return Q(I,0.,(1./3),0.,q)/Q(I,0.,(1./3),0.,1.)

    class Message():

        def __init__(self,N,I,nplot=1000):

            self.N = N
            self.I = I
            self.M = D(I,self.N) 
            self.nplot = nplot

        def error(self):
            """Returns the L^{2} error"""
            E = 0.
            I0 = integrate.quad(lambda x: C(self.I,x)**2,0.,1.)[0]
            for n in range(0,self.N):
                f = lambda x: (n/(1.*self.N)-C(self.I,x))**2
                E = E + integrate.quad(f,self.M[n],self.M[n+1])[0]
            return  np.sqrt(E/I0)

    def mplot(M,nplot=1000):
        """Plots a list of message functions"""
        I = M[0].I
        q = np.linspace(0.,1.,nplot)
        cplt = np.zeros(nplot)
        for n in range(0,nplot):
            cplt[n] = C(I,q[n])
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
        """Writes a list of message function to disk"""
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
                R = [q[n],I(q[n]),C(I,q[n])]
                for m in M:
                    R = R + [np.digitize(q[n],m.M)/(1.*m.N)] 
                writer.writerow(R)
                

    def eplot(I,nplot=4):
        """Plots the error against the number of messages"""
        E = np.zeros(nplot)
        for n in range(0,nplot):
            E[n] = Message(10**(n+1),I).error()
        plt.plot(np.r_[1:(nplot+1)],np.log10(E))
        plt.xlabel('$\\log_{10}$(Number of Messages)')
        plt.ylabel('$\\log_{10}$(Relative $L^{2}$ Error)')
        plt.grid()
        plt.show()

def _fixed_quad(*args,**kwargs):
    """Wrapper for fixed_quad that returns the integral (without the error)"""
    return fixed_quad(*args,**kwargs)[0]

# -----------------------------------------------------------------------------
# functions
def identity_cost(func,e_bar,tol=1.e-10):
    """
    Cost of identity message function
    
    Parameters
    ----------
    func : callable
        PDF of error. identity_cost assumes func has support [-e_bar,e_bar]
    e_bar : float
        Support parameter for the error distribution (see above).

    Notes
    -----
    The action function kinks at 1-e_bar and 1+e_bar. Therefore, identity_cost 
    breaks the integration problem into three parts. The action function is 
    computed using scipy.integrate._fixed_quad, which caches weights/knots, so
    should be fast. 

    There are 9 anonymous functions. Guido would not approve.
    """

    if np.abs(_fixed_quad(func,-e_bar,e_bar)-1.)>tol:
        raise Exception('PDF of error does not integrate to one')

    # quasi-expectation
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

    # integration over [+e_bar,1-e_bar]
    a_md = lambda m_til: A(m_til,m_til-e_bar,m_til+e_bar)
    f_md = lambda q, m_til: func(m_til-q)*(q-a_md(m_til))**2.
    z_md = dblquad(f_md,e_bar,1.-e_bar,gfun,hfun)

    # integration over [1-e_bar,1+e_bar]
    a_up = lambda m_til: A(m_til,m_til-e_bar,1.)
    f_up = lambda q, m_til: func(m_til-q)*(q-a_up(m_til))**2.
    z_up = dblquad(f_up,1.-e_bar,1.+e_bar,gfun,1.)

    # integration over [-e_bar,+e_bar]
    a_dn = lambda m_til: A(m_til,0.,m_til+e_bar)
    f_dn = lambda q, m_til: func(m_til-q)*(q-a_dn(m_til))**2.
    z_dn = dblquad(f_dn,-e_bar,e_bar,0.,hfun)

    return z_up[0] + z_md[0] + z_dn[0]

def discrete_cost(func,N,tol=1.e-10):
    """
    Cost of discrete message function
    
    Parameters
    ----------
    func : callable
        PDF of error distribution. identity_cost assumes func has support
        [-e_bar,e_bar] where e_bar = 1/2N 
    N : int
        There are N+1 messages
    """

    e_bar = 1./(2.*N) 
    
    if np.abs(_fixed_quad(func,-e_bar,e_bar)-1.)>tol:
        raise Exception('PDF of error does not integrate to one')

    _exp = _fixed_quad(lambda e : func(e)*e,-e_bar,e_bar)
    _var = _fixed_quad(lambda e : func(e)*(e-_exp)**2.,-e_bar,e_bar)
    return (N+1.)*_var
