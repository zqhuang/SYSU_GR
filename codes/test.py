import sympy as sym
x, y = sym.symbols('x, y')
f = x**2
y = f+x
print(sym.diff(y, x))
