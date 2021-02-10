## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
#print(sym.integrate(sym.log(x),x))
print(sym.integrate(sym.sinh(x)**2,x))
#print(sym.integrate(sym.exp(-x**2)*sym.erf(x),x))
#print(sym.integrate(sym.exp(-x**2),(x, -sym.oo, sym.oo)))
