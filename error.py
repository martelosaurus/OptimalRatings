import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class error:
    """
    Error RV (Beta(a,b) shifted from [0,1] to [-e_bar,+e_bar]) 
    RV has generic name "X" in documentation

    TODO: subclass scipy.stats.rv_continuous (seems complicated a.f.) 
    """

    def __init__(self,a,b,e_bar):
        """
        Parameters
        ----------
        a : float 
            Beta parameter (a>0)
        b : float 
            Beta parameter (b>0)
        e_bar : float
            Support is [-e_bar,+e_bar]

        Functions and Parameters
        ------------------------
        __init__ constructs the following functions and parameters
        x : function
            Maps x\in[0,1] to e\in[-\\bar{e},+\\bar{e}]
        f : function 
            Non-normalized Beta(a,b) PDF 
        n : float
            Normalization

        TODO: assumes default options are fine
        """
        self.a = a
        self.b = b
        self.e_bar = e_bar

        x = np.vectorize(lambda e : .5*(1.+e/e_bar))
        self.f = np.vectorize(lambda e : x(e)**(a-1.)*(1.-x(e))**(b-1.))
        self.n = quad(self.f,-e_bar,e_bar) 

    @np.vectorize
    def pdf(self,x):
        """(p)robability (d)ensity (f)unction"""
        return self.f(x)/self.n

    def expect(self,func,lb,ub):
        """
        Expectation of X given the condition lb < X < ub

        Parameters
        ----------
        func: callable
            Function for which integral is calculated
        lb : float
            Lower limit of condtion 
        ub : float
            Upper limit of condtion 

        TODO: assumes default options are fine
        """
        If = quad(lambda x: func(x)*self.f(x),lb,ub)[0]
        I0 = quad(lambda x: self.f(x),lb,ub)[0]
        return If/I0

    def pdf_plot(self):
        """
        Plots the PDF of the distribution

        Returns
        -------
        A plt object
        """
        return 


