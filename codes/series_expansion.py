## example of differentiation
import sympy as sym
x = sym.Symbol('x')
#print(sym.series(sym.cos(x),x))
#print(sym.series(sym.cos(x),x, 0, 12))
print(sym.series(sym.sin(sym.sinh(x))-sym.sinh(sym.sin(x)),x,0,8))
