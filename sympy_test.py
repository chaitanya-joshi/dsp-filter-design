from numpy import zeros, ceil, sin
from sympy import symbols, expand, Function
from sympy.plotting import plot
# N/2from sympy import *
s, x = symbols('s x')
expr = x**2
expr = expr.subs(x,1j*x)
plot(expr, (x, -5, 5))
# transform = (s**2 + 1)/s
# bilinear = (1.0 - x)/(1.0 + x)
# expr = 1.0/(1.0 + s**2)
# print expr
# expr2 = expr.subs(s, transform)
# print expr2
# expr3 = expr2.subs(s, bilinear)
# print expr3
# expr2 = (1 + x**3)/((1+x)*(1+2*x))
# expr3 = expand(x*expr2)
# factors = 1.0
# factors = factors*expr2
# sq =  expand(pow((1.0+x),3))
# f = Function('f')(x)
# f = x + 2*x**2
# # k = pow(y,-1)
# # print expr2.subs(x,k)
# print cancel(expr2)
#
# a = [1,2,3]
# print a[-1]
# N = 10.0
# k = ceil(N/2)
# ax = [1.0,2.0] # Array of size [(order of ax) + 1] with ax[0] being coeff of x^n
# ax = list(reversed(ax))
# ax = [i/2.0 for i in ax]
# print ax