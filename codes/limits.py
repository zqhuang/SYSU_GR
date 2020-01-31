## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
print(sym.limit((sym.log(1+x)-x)/x**2, x, 0))
