import numpy as np

M = 4
N = 2
K = M*N
e_bar = .5/N
d_bar = .5/K

def k_up(i):
    return min(i+M,K)+1

def k_dn(i):
    return max(0,i+1)

# A - fun
def A(i,j):
    return _A[i+M,j]

# A
a_str = []
a_file = open("a_file.tex","w")
for i in range(-M,K):
    row = []
    for j in range(0,K+1):
        if j-M<=i<j:
            row.append("\\alpha_{"+str(i)+","+str(j)+"}")
        else:
            row.append("")
    a_str.append(" & ".join(row))
a_str = "\\\\\n".join(a_str)
a_str_up = "\\begin{array}{" + (K+2)*"c" + "}"
a_str_dn = "\\end{array}"
a_str = "\n".join([a_str_up,a_str,a_str_dn])
a_file.write(a_str)
a_file.close()

# B
b_str = []
b_file = open("b_file.tex","w")
for i in range(-M,K):
    row = []
    for j in range(1,K+1):
        if j-1-M == i:
            #B[i+M,j-1] = A(j-1-M,j-1)
            row.append("\\alpha_{"+str(j-1-M)+","+str(j-1)+"}")
        elif j-1-M<i<j-1:
            #B[i+M,j-1] = A(i,j-1)-A(i,j)
            row.append("\\alpha_{"+str(i)+","+str(j-1)+"}-\\alpha_{"+str(i)+","+str(j)+"}")
        elif i == j-1:
            #B[i+M,j-1] = -A(j-1,j)
            row.append("-\\alpha_{"+str(j-1)+","+str(j)+"}")
        else:
            #pass
            row.append("")
    b_str.append(" & ".join(row))
b_str = "\\\\\n".join(b_str)
b_str_up = "\\begin{array}{" + K*"c" + "}"
b_str_dn = "\\end{array}"
b_str = "\n".join([b_str_up,b_str,b_str_dn])
b_file.write(b_str)
b_file.close()

# C 
c_str = []
c_file = open("c_file.tex","w")
for i in range(-M,K):
    row = []
    for j in range(0,K+2):
        if k_dn(i) == j:
            #C[i+M,j] = -A(i,k_dn(i))
            row.append("-\\alpha_{"+str(i)+","+str(k_dn(i))+"}")
        elif k_dn(i)<j<k_up(i):
            #C[i+M,j] = A(i,j-1)-A(i,j)
            row.append("\\alpha_{"+str(i)+","+str(j-1)+"}-\\alpha_{"+str(i)+","+str(j)+"}")
        elif j == k_up(i):
            #C[i+M,j] = A(i,k_up(i)-1)
            row.append("\\alpha_{"+str(i)+","+str(k_up(i)-1)+"}")
        else:
            #pass
            row.append("")
    c_str.append(" & ".join(row))
c_str = "\\\\\n".join(c_str)
c_str_up = "\\begin{array}{" + (K+2)*"c" + "}"
c_str_dn = "\\end{array}"
c_str = "\n".join([c_str_up,c_str,c_str_dn])
c_file.write(c_str)
c_file.close()
