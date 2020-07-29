class Identity(Message):

    def __init__(self):

    def cost(self.e_bar,tol=1.e-10):
        """
        Cost of identity message function with constant importance function

        'identity_cost' stands alone because it needs to run faster than 
        Smooth.cost for many applications.  
        
        Parameters
        ----------
        func : callable
            PDF of error. identity_cost assumes func has support [-e_bar,e_bar]
        self.e_bar : float
            Support parameter for the error distribution (see above).
        tol : float
            Error PDF should integrate to unity (+/- tol)

        Notes
        -----
        The action function kinks at 1-self.e_bar and 1+self.e_bar. Therefore, 
        identity_cost breaks the integration problem into three parts. The 
        action function is computed using scipy.integrate.fixed_quad, which 
        caches weights/knots, so should be fast. 

        Note if the error is uniformly distributed, then 

            cost[identity] = (1./3.)*(e_bar**2.)*(1.-e_bar)

        There are 9 anonymous functions. Guido would not approve.
        """

        if np.abs(_fixed_quad(func,-self.e_bar,self.e_bar)-1.)>tol:
            raise Exception('PDF of error does not integrate to one')

        # quasi-expectation
        @np.vectorize
        def A(m_til,a,b):
            if a == b:
                return a
            else:
                I0 = _fixed_quad(lambda x: self.func(m_til-x)*x**0.,a,b)
                I1 = _fixed_quad(lambda x: self.func(m_til-x)*x**1.,a,b)
                return I1/I0

        # constant m_tilde limits of integration (see dblquad doc)
        afun = lambda _ : 0.
        bfun = lambda _ : 1.

        # non-constant limits 
        gfun = lambda m_til : m_til-self.e_bar # lower
        hfun = lambda m_til : m_til+self.e_bar # upper 

        # integration over [1-self.e_bar,1+self.e_bar]
        a_up = lambda m_til: A(m_til,m_til-self.e_bar,1.)
        f_up = lambda q, m_til: self.func(m_til-q)*(q-a_up(m_til))**2.
        z_up = fixed_dblquad_tri(f_up,1.-self.e_bar,1.+self.e_bar,gfun,bfun)

        # integration over [+self.e_bar,1-self.e_bar]
        a_md = lambda m_til: A(m_til,m_til-self.e_bar,m_til+self.e_bar)
        f_md = lambda q, m_til: self.func(m_til-q)*(q-a_md(m_til))**2.
        z_md = fixed_dblquad_tri(f_md,self.e_bar,1.-self.e_bar,gfun,hfun)

        # integration over [-self.e_bar,+self.e_bar]
        a_dn = lambda m_til: A(m_til,0.,m_til+self.e_bar)
        f_dn = lambda q, m_til: self.func(m_til-q)*(q-a_dn(m_til))**2.
        z_dn = fixed_dblquad_tri(f_dn,-self.e_bar,self.e_bar,afun,hfun)

        return z_up + z_md + z_dn
