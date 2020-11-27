## example of differentiation
import sympy as sym
x, y = sym.symbols('x,y')
A = sym.Matrix([[1,x], [y,1]])
print(A**(-1))
