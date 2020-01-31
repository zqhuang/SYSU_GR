## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
f = x ** 4 - 3 * x ** 2 + 1
print(sym.factor(f))
print(sym.factor(f, modulus=5))
