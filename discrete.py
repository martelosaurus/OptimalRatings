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

    def cost(self):
        """
        Cost of discrete message function 
        
        Parameters
        ----------
        N : int
            There are N+1 messages

        """
        return 1. 
