## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y, z = sym.symbols('y, z')
print(sym.diff(sym.log(x)))
print(sym.diff(x*sym.exp(x+y), x))
