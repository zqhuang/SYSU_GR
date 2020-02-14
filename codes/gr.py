from gravipy.tensorial import * # import GraviPy package
from sympy import init_printing
import inspect
init_printing()

##-----------------------------------------------------------
#Each component of any tensor object, can be computed by calling the appropriate instance of the GeneralTensor subclass with indices as arguments. The covariant indices take positive integer values (1, 2, ..., dim). The contravariant indices take negative values (-dim, ..., -2, -1).

t, r, theta, phi, M = symbols('t, r, \\theta, \phi, M')
# create a coordinate four-vector object instantiating  the Coordinates class
x = Coordinates('\chi', [t, r, theta, phi])
# define a matrix of a metric tensor components
Metric = diag(-(1-2*M/r), 1/(1-2*M/r), r**2, r**2*sin(theta)**2)  
# create a metric tensor object instantiating the MetricTensor class
g = MetricTensor('g', x, Metric)


Gam = Christoffel('Gam', g)

print(Gam(-1, 2, 1))
