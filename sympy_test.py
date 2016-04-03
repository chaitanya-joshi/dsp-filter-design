from numpy import zeros
# from sympy import symbols, expand, Function
from sympy import *
x, y = symbols('x y')
expr = 1.0 + 2*y
expr2 = x + 2.12345*x**2
expr3 = expand(x*expr2)
factors = 1.0
factors = factors*expr2
sq =  expand(pow((1.0+x),3))
f = Function('f')(x)
f = x + 2*x**2
k = pow(y,-1)
print expr2.subs(x,k)