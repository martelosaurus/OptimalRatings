def f(n,N):
	if n == N:
		return N
	else:
		return n*(1+f(n-1,N))
