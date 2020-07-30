# --------------------------------------------------------------------------- #
import numpy as np
from scipy import linalg as la
from utils import beta

M = 2
N = 2
K = M*N
e_bar = .5/N
d_bar = .5/K

def k_up(i):
    return min(i+M,K)+1

def k_dn(i):
    return max(0,i+1)

# --------------------------------------------------------------------------- #
if False:
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

    if True:
        for i, j in zip(z1, z2):
            print(i + j)

# --------------------------------------------------------------------------- #
# A
#A = np.random.randint(1,10,(M+K,K+2));
#c = np.random.randint(1,10,M+K);
#r = np.random.randint(1,10,K+2);
if False:
    c = np.ones((1,M+K))
    r = np.ones((1,K+2))
    _A = la.hankel(c,r)
_A = np.zeros((M+K,K+2))
for i in range(-M,K):
    for j in range(0,K+2):
        if j-M<=i<j:
            _A[i+M,j] = 1./M
def A(i,j):
    return _A[i+M,j]

# C 
if True:
    C = np.zeros((M+K,K+2))
    for j in range(0,K+2):
        for i in range(-M,K):
            if k_dn(i) == j:
                C[i+M,j] = -A(i,k_dn(i))
            elif k_dn(i)<j<k_up(i):
                C[i+M,j] = A(i,j-1)-A(i,j)
            elif j == k_up(i):
                C[i+M,j] = A(i,k_up(i)-1)
            else:
                pass

# B
B = np.zeros((M+K,K))
for i in range(-M,K):
    for j in range(1,K+1):
        if j-1-M == i:
            B[i+M,j-1] = A(j-1-M,j-1)
        elif j-1-M<i<j-1:
            B[i+M,j-1] = A(i,j-1)-A(i,j)
        elif i == j-1:
            B[i+M,j-1] = -A(j-1,j)
        else:
            pass
