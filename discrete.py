import utils 

@np.vectorize
def _alpha(i,j,e_bar):
    return fixed_quad(lambda e : _beta(e,a,b,e_bar),-e_bar,e_bar)[0]

class Discrete(Message):

    def __init__(self,N,M):

        # number of knots
        self.M = M
        self.N = N

        # knot spacing
        self.e_bar = .5/N
        self.d_bar = .5/K

        # knots
        self.y

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


    def cost(self):
        """
        Cost of discrete message function 
        
        Parameters
        ----------
        N : int
            There are N+1 messages

        """
        return 1. 
