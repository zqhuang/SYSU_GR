## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
print(sym.diff(sym.sin(x**2*y+y),x))
print(sym.diff(sym.tan(x)))
