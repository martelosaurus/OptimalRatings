from discrete import Message
from matplotlib import pyplot as plt
import pandas as pd

X = pd.DataFrame()
Y = pd.DataFrame()
B = [-.2,-.1,0.,1.,2.,4.]
while B:
    n = len(B)
    _b = B.pop()
    m = Message(M=100,N=4,b=_b)
    xx, mm = m.plot_msg("msg" + str(n) + ".pdf",title=False)
    ee, yy = m.plot_err("err" + str(n) + ".pdf",2.,6.,title=False)
    X["x" + str(n)] = xx
    X["m" + str(n)] = mm
    Y["e" + str(n)] = ee
    Y["y" + str(n)] = yy
X.to_csv("msg.csv")
Y.to_csv("err.csv")
