import sympy as sym

a, b, x = sym.symbols('a b x') 
a = 1-b/2

sol = sym.idiff(2*x*(6*a+3*b*x)*(6*a+3*b*(1+x))-(3*a*x+2*b*x**2)*(6*a+3*b*(1+x))-(3*a*(1+x)-2*b*(x**2+x+1))*(6*a+3*b*x),x,b)

