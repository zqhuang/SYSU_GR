## example of differentiation
import sympy as sym
import gravipy as gr
x, y = sym.symbols('x, y')
f, g = sym.symbols('f, g', cls = sym.Function)
print(sym.dsolve( f(x).diff(x, x) + f(x), f(x) ))
print(sym.dsolve(sym.sin(x) * sym.cos(f(x)) + sym.cos(x) * sym.sin(f(x)) * f(x).diff(x), f(x), hint='separable'))
print(sym.dsolve( x * f(x).diff(x) + f(x) - f(x) ** 2, f(x)))
print(sym.dsolve( x * f(x).diff(x) + f(x) - f(x) ** 2, f(x), hint='Bernoulli'))
