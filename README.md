# Ratings and Cooperative Information Transmission

Forthcoming at *Management Science*

We look at two dimensions: the importance function and the error distribution. In one set of analyses, we assume uniform importance and non-uniform error. In the other, we assume non-uniform importance and uniform error.

## Non-Uniform error, Uniform Importance

Discrete messages

```python
from discrete import Message
from matplotlib import pyplot as plt

B = [-.2,-.1,0.,1.,2.,4.]
while B:
    n = len(B)
    _b = B.pop()
    m = Message(M=100,N=4,b=_b)
    m.plot_msg("msg" + str(n) + ".pdf",title=False)
    m.plot_err("err" + str(n) + ".pdf",2.,6.,title=False)
```

## Uniform Error, Non-Uniform Importance

```python
from optrat import Message, mplot
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
        M = M + [Message(n,I[i])]
   	mplot(M)
```
