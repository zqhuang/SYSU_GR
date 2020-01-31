## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
print(sym.series(sym.cos(x),x))
print(sym.series(sym.cos(x),x, 0, 12))
print(sym.series(2/x-2/sym.cosh(x)/x,x, 0, 7))
