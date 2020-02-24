import sympy as sym

t = sym.Symbol('t')

vec = [ t**2, t ]

def d_dt(vec):
    global t
    return [ sym.simplify(sym.diff(x, t)) for x in vec ]

print d_dt(vec)
