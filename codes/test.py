import numpy as np
dist = 3.
q_arr = [1.]
z_arr = [dist]

def add_point():
    global q_arr, z_arr, dist
    n = len(q_arr)
    a = dist + z_arr[n-1]
    q_arr.append(q_arr[n-1]/a)
    z_arr.append(dist-1/a)

def potential(x, y, z):
    global q_arr, z_arr, dist
    p = 0.
    for i in range(len(q_arr)):
        p = p + q_arr[i]*(1./np.sqrt((z_arr[i]-z)**2+x**2+y**2)-1./np.sqrt(x**2+y**2+(z_arr[i]+z)**2))
    p = p / sum(q_arr)
    return p
        
for i in range(12):
    add_point()
    

print(potential(0., 0.,dist/2.)/potential(0., 0., dist+1.))
print((2./dist - 1./(1.5*dist))/(1.-1./2./dist))
