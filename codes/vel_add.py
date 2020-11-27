## sample script of Lorentz Transformation
## use natural units
import sympy as sym

def add_velocity_vector(u, v, theta):
    return [ (u[0]+v)/(1+u[0]*v), u[1]*sym.sqrt(1-v**2)/(1+u[0]*v), u[2]*sym.sqrt(1-v**2)/(1+u[0]*v) ]

def add_velocity(u, v, theta):
    return sym.sqrt(u**2+v**2+2*u*v*sym.cos(theta)-(u*v*sym.sin(theta))**2)/(1+u*v*sym.cos(theta))


    
#check the formulas
u, v, theta = sym.symbols('u, v, theta')
uvec = [ u*sym.cos(theta), u*sym.sin(theta), 0 ]
svec = add_velocity_vector(uvec, v, theta)
print(sym.simplify(svec[0]**2+svec[1]**2+svec[2]**2-add_velocity(u, v, theta)**2))
print(add_velocity(0.2, 0.2, sym.pi/3).evalf(5))
