## sample script for frenet calculus
#----------by Zhiqi Huang--for the course ''General Relativity''-------------

import sympy as sym

t = sym.symbols('t') #curve parameter

####################for a general test######################################
####This somehow does not work; guess because sympy sucks on abstract derivation
#############################################################################
#x, y, z = sym.symbols('x, y, z', cls = sym.Function)
#r = [ x(t), y(t), z(t) ]


####################testing a specific model #####################
#This works perfectly fine for simple models #####################
#For more complicated models you may want to switch off check_s###
check_s = False #calculate everything by converting t to s (length parameter)
lam = sym.Symbol('lambda', positive=True)
r = [ sym.cos(t), sym.sin(t), t ]
###################################################################

rp = [sym.diff(r[0],t), sym.diff(r[1],t), sym.diff(r[2],t)]  ## this is d r / dt
def magnitude(vec):
    return sym.sqrt(sum([x**2 for x in vec]))

dsdt = magnitude(rp)

def normalize(vec):
    norm = magnitude(vec)
    return [ sym.simplify(x/norm) for x in vec ]

def linear_add(a, b, ca, cb):
    return [ sym.simplify(ca*a[0]+cb*b[0]), sym.simplify(ca*a[1]+cb*b[1]), sym.simplify(ca*a[2]+cb*b[2]) ]

def prod_outer(a, b):
    return [ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] ]

def prod_inner(a, b):
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def dvecdt(vec):
    global t
    return [ sym.diff(x, t) for x in vec ]

def iszero(vec):
    return vec[0].equals(0) and vec[1].equals(0) and vec[2].equals(0)

def tangent_vector():
    global rp
    return normalize(rp)

##=================first I derive everything with s-based formulas=================
def dvecds(vec):
    global dsdt
    dvdt = dvecdt(vec)
    return [ vp/dsdt for vp in dvdt ]

def normal_vector_s():
    return normalize( dvecds(tangent_vector()) )

def curvature_s():
    return sym.simplify(magnitude( dvecds(tangent_vector()) ))

def tortion_s():
    rs = tangent_vector()
    rss = dvecds(rs)
    rsss = dvecds(rss)
    vol = prod_inner(prod_outer(rs, rss), rsss)
    kappa = magnitude(rss)
    return sym.simplify(vol/kappa**2)

def binormal_vector_s():
    rs = tangent_vector()
    rss = dvecds(rs)
    return normalize( prod_outer(rs, rss) )

##===============================now I enter the t-based formulas===============================
def normal_vector_t():
    global rp
    rpp = dvecdt(rp)
    mag = magnitude(prod_outer(rp, rpp))
    return linear_add(rpp, rp, dsdt/mag, -prod_inner(rp, rpp)/dsdt/mag)

def binormal_vector_t():
    global rp
    rpp = dvecdt(rp)
    return normalize(prod_outer(rp, rpp))

def curvature_t():
    global rp
    rpp = dvecdt(rp)
    p = prod_outer(rp, rpp)
    return sym.simplify(magnitude(p)/dsdt**3)
    
def tortion_t():
    global rp
    rpp = dvecdt(rp)
    rppp = dvecdt(rpp)
    p = prod_outer(rp, rpp)
    return sym.simplify(prod_inner(prod_outer(rp, rpp), rppp)/magnitude(p)**2)

print("-----------tangent vector--------------")
print(tangent_vector())
print("-----------normal vector---------------")
print(normal_vector_t())
if(check_s):
    print(normal_vector_s())
print("-----------binormal vector-------------")
print(binormal_vector_t())
if(check_s):
    print("check by using s variable:")    
    print(binormal_vector_s())
print("-----------curvature-------------------")
print(curvature_t())
if(check_s):
    print("check by using s variable:")
    display(curvature_s())
print("-----------tortion---------------------")
print(tortion_t())
if(check_s):
    print("check by using s variable:")    
    display(tortion_s())


