import utils 

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

    def cost(self,N):
        """
        Cost of discrete message function 
        
        Parameters
        ----------
        N : int
            There are N+1 messages

        """
        return 1. 
