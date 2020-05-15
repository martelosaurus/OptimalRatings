"""
This is a sandbox script for the discrete message module
"""
import numpy as np
import dm as dm

# importance function(s)
I = {
        'i1' : lambda x: x**3., 
        'i2' : lambda x: x**(-1.5),
        'i3' : lambda x: (6.*(x-.5)**2.+.5)**3.
    }

# for each importance function, plot discrete messages of size 5, 20, 100 
for i in I:
    M = []
    for n in [5,20,100]:
        M = M + [dm.Message(n,I[i])]
    dm.mplot(M)
    
