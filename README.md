# Ratings and Cooperative Information Transmission

Forthcoming at *Management Science*

The sender and receiver agree on a message function `m:[0,1]->[0,1]`. The sender privately observes `q`, which is uniformly distributed on `[0,1]`. She sends the receiver the message `m(q)`. The receiver receives the message `m_tilde = m(q)+e`, where `e` is distributed on `[-e_bar,e_bar]` according to the PDF `f`. She then takes an action `A(m_tilde)`. The sender and receiver incur the cost `((q-A(m_tilde))^2)*I(q)` where `I:[0,1]->[0,1]` is the importance function. 

We look at two dimensions: the importance function and the error distribution. In one set of analyses, we assume uniform importance and non-uniform error. In the other, we assume non-uniform importance and uniform error. 

## Non-Uniform Error, Uniform Importance

`I(q)=1` and `f` is the PDF of a random variable that is distributed according to a Beta distribution on `[-e_bar,e_bar]`. One of the two parameters is pinned-down by the assumption that `f` is symmetric about zero. The other parameter is free and denoted by `b`.


```python
from constant_I import Message

B = [-.2,-.1,0.,1.,2.,4.]
while B:
    n = len(B)
    _b = B.pop()
    m = Message(M=100,N=4,b=_b)
    m.plot_msg("msg" + str(n) + ".pdf",title=False)
    m.plot_err("err" + str(n) + ".pdf",2.,6.,title=False)
```

## Uniform Error, Non-Uniform Importance

`f(e)=1/(2*e_bar)` is constant and `I` is specified below.

```python
from constant_f import Message

# importance function(s)
I = {
    'i1' : lambda x: x**3., 
    'i2' : lambda x: x**(-1.5),
    'i3' : lambda x: (6.*(x-.5)**2.+.5)**3.
    }

# for each importance function, plot discrete messages of size 5, 20, 100 
for i in I:
    for n in [5,20,100]:
        m = Message(n,I[i])
        m.plot_msg("msg" + i + str(n) + ".pdf",title=False)
```
