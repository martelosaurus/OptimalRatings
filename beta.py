import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta

# -----------------------------------------------------------------------------
# parameters
n_plot = 100 
a_max = 2.
b_max = 2.

# -----------------------------------------------------------------------------
# plotting
A, B = np.meshgrid(np.linspace(0.,a_max,n_plot),np.linspace(0.,b_max,n_plot))




