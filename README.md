# Ratings and Cooperative Information Transmission
Numerical solutions for Optimal Ratings

```python
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
```

Discrete messages

```python
from discrete import Message

m = Message(M=100,N=4,b=2)
m.plot_msg("../new_tex/Model_Figures/msg" + str(n) + ".pdf",title=False)
m.plot_err("../new_tex/Model_Figures/err" + str(n) + ".pdf",2.,6.,title=False)
```
