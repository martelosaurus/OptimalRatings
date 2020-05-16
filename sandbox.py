from matplotlib import pyplot as plt
import numpy as np

n_plot = 100
x_plot = np.linspace(0.,1.,n_plot)
y_plot = np.linspace(0.,1.,n_plot)
X, Y = np.meshgrid(x_plot,y_plot)

@np.vectorize
def f(x,y):
    if x>y:
        return np.nan
    else:
        return x**2.+y**2.

plt.contourf(X,Y,f(X,Y))
plt.show()
