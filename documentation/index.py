# --------------------------------------------------------------------------- #
import numpy as np

r = 2
M = 4
N = 5
K = M*N

def ran(a,b):
    return range(a,b+1)

def k_up(i):
    return min(i+M,K)

def k_dn(i):
    return max(0,i+1)

Z = np.random.rand(M+K-1,K+1);
Z = Z.round(r)

# method 2
z2 = []
for j in ran(0,K):
    for i in ran(j-M,j-1):
        z2.append(str([i,j]))
z2.sort()

# method 1
z1 = []
for i in ran(-M,K-1):
    for j in ran(k_dn(i),k_up(i)):
        z1.append(str([i,j]))
z1.sort()

# difference
print(set(z1) - set(z2))
print(set(z2) - set(z1))

if False:
    for i, j in zip(z1, z2):
        print(i + j)


