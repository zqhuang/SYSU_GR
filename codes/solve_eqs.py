## example of differentiation
import sympy as sym
x = sym.Symbol('x')
y = sym.Symbol('y')
print(sym.solveset(x**3+1, x))
print(sym.solveset(sym.exp(x)+1, x))
solution = sym.solve( (x + 5 * y - 2, -3 * x + 6 * y - 15), (x,y))
print(solution)
print(solution[x])
print(solution[y])
