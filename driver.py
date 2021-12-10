# ------------------------------------------------------------------------------
# FIGURE 6
from constant_I import Message

B = [-.2,-.1,0.,1.,2.,4.]
while B:
    n = len(B)
    _b = B.pop()
    m = Message(M=100,N=4,b=_b)
    m.plot_msg("msg" + str(n) + ".pdf",title=False)
    m.plot_err("err" + str(n) + ".pdf",2.,6.,title=False)

# ------------------------------------------------------------------------------
# FIGURES 8-10
from constant_f import Message

# importance function(s)
I = {
    'i1' : lambda x: x**3., 
    'i2' : lambda x: x**(-1.5),
    'i3' : lambda x: (6.*(x-.5)**2.+.5)**3.
    }

# for each importance function, plot discrete messages of size 5 and 20
for i in I:
    for n in [5,20]:
        m = Message(n,I[i])
        m.plot_msg("msg" + i + str(n) + ".pdf",title=False)

# ------------------------------------------------------------------------------
# FIGURE 11: DJIA

# Data Note: Dow Jones Industrial Average, 1885-02-16 through 20210-09-03, retrieved from MeasuringWorth; https://www.measuringworth.com/datasets/DJA/index.php, September 7, 2021. 
from djia import _djia
_djia()

# ------------------------------------------------------------------------------
# FIGURES 12-13: WFA

# Data Note: Courtesy of Joseph Zechner.
from wfa import _wfa
_wfa()

